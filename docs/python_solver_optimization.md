# Python Solver Optimization Review

## Executive Summary

| Component | Current Time | Target | Strategy |
|-----------|--------------|--------|----------|
| C++ UMAT call | 4.8 µs | 4.8 µs | Already optimized (zero-copy views) |
| np.linalg.solve | 3.3 µs | 3.3 µs | NumPy overhead unavoidable |
| Python overhead | 6.4 µs | 2-3 µs | See options below |
| **Total** | **14.5 µs** | **10-11 µs** | 1.3x speedup possible |

**Performance vs C++:** Currently 1.2x ratio (14.5 µs Python vs ~12 µs C++)

The goal is to reduce Python overhead from ~6 µs to 2-3 µs while maintaining the flexibility of the Python API.

---

## 1. Current Bottleneck Analysis

### Detailed Profiling Results

**Measured using `timeit` with 10,000+ iterations:**

| Operation | Time (µs) | % of Total | Notes |
|-----------|-----------|------------|-------|
| `scc.umat_inplace()` | 4.77 | 33% | C++ UMAT + binding overhead |
| `np.linalg.solve(6x6)` | 3.29 | 23% | Newton-Raphson linear solve |
| `HistoryPoint.from_state()` | 1.31 | 9% | History storage |
| `reshape` views (×2) | 0.56 | 4% | Array reshaping for binding |
| Python solver loop | 4.57 | 31% | Control flow, attribute access |
| **Total** | **14.5** | 100% | |

**Python vs C++ breakdown:**
- Pure computation: ~8 µs (55%)
- Python overhead: ~6.5 µs (45%)

### Time Breakdown (per increment)

```
Total: 14.5 µs
├── C++ umat_inplace:     4.8 µs (33%)  - Already optimized with carma views
├── np.linalg.solve:      3.3 µs (23%)  - Mixed boundary conditions
├── History storage:      1.3 µs (9%)   - HistoryPoint.from_state()
├── Python solver loop:   4.6 µs (32%)  - Main optimization target
│   ├── Function calls:   ~1.5 µs
│   ├── Attribute access: ~1.2 µs
│   ├── Dataclass ops:    ~1.0 µs
│   └── Control flow:     ~0.9 µs
└── Array reshaping:      0.5 µs (3%)   - reshape for C++ binding
```

### Profiling the Solver Loop

```python
# Key operations in _solve_increment (strain-controlled case):
def _solve_increment(...):
    # 1. Apply strain increment
    sv.DT = Dtinc * DT_target          # Attribute access
    sv.DR.fill(0.0)                     # NumPy method call
    np.fill_diagonal(sv.DR, 1.0)        # NumPy function call
    np.copyto(sv.DEtot, DEtot_target)   # NumPy function call
    sv.DEtot *= Dtinc                   # NumPy in-place op

    # 2. Call UMAT
    self._call_umat(block, sv, ...)     # Method call + C++ binding

    # 3. Return
    return True
```

Each Python function call has ~50-100 ns overhead. With 20+ operations per increment, this adds up.

---

## 2. Binding Layer Analysis: Armadillo+Carma vs Eigen

### 2.1 Current Architecture

The current simcoon Python binding uses:
- **Armadillo** (C++ linear algebra library) for all matrix operations
- **Carma** (Armadillo ↔ NumPy bridge) for zero-copy data transfer
- **pybind11** for Python binding generation

**Binding call path:**
```
Python (NumPy array)
    ↓ (zero-copy via carma::arr_to_mat_view)
pybind11 dispatch
    ↓ (function lookup by name)
Armadillo wrapper
    ↓ (matrix operations)
UMAT implementation
    ↓ (stress/tangent computation)
Return via carma view (zero-copy)
```

### 2.2 Carma Zero-Copy Analysis

Carma already provides **optimal zero-copy transfer**:

```cpp
// Current implementation in run_umat.cpp
void umat_inplace(const string& umat_name,
                  py::array_t<double>& etot,  // NumPy array
                  ...) {
    // Zero-copy view (no data copy!)
    auto etot_mat = carma::arr_to_mat_view(etot);

    // Direct operation on NumPy memory
    umat_dispatcher(umat_name, etot_mat, ...);

    // Changes visible in Python immediately
}
```

**Carma overhead measured: ~0.3-0.5 µs per array view creation**

This is already near-optimal. Further reduction would require:
- Removing pybind11 dispatch overhead (~0.1 µs)
- Removing UMAT name lookup overhead (~0.05 µs)

### 2.3 Eigen Migration Assessment

**Would migrating to Eigen+pybind11 (without Carma) help?**

| Factor | Armadillo+Carma | Eigen+pybind11 | Winner |
|--------|-----------------|----------------|--------|
| Zero-copy to NumPy | Yes (carma views) | Yes (Eigen::Map) | Tie |
| Binding overhead | ~0.4 µs | ~0.3 µs | Eigen (marginal) |
| Code readability | MATLAB-like | Template-heavy | Armadillo |
| Existing codebase | 226 files, 1267 usages | 0 | Armadillo |

**Migration effort estimate:**
- 226 C++ files use Armadillo
- 1,267 `arma::mat` occurrences
- 919 `arma::vec` occurrences
- 160+ files in Continuum_mechanics/
- 41 UMAT models

**Estimated effort: 1,900 - 2,950 hours (50-75 person-weeks)**

**Recommendation: Do NOT migrate to Eigen.**

The binding overhead difference is marginal (~0.1 µs), while migration cost is enormous. The real bottleneck is the Python loop, not the C++ binding.

### 2.4 Alternative: Nanobind

If binding overhead becomes critical, consider **nanobind** (pybind11 successor):
- 2-3x faster dispatch than pybind11
- Smaller binary size
- Drop-in replacement for most pybind11 code

**Estimated binding overhead reduction: 0.1-0.2 µs per call**

However, this requires:
- Updating all binding code
- Testing compatibility with carma
- May not work with all carma features

**Recommendation:** Keep as future option if Python loop is optimized first.

---

## 3. Optimization Strategies

### Strategy A: Cython Solver Loop (Recommended)

**Effort: Medium | Impact: High (2-3x speedup)**

Convert the hot path (`_solve_increment`, `_call_umat`) to Cython while keeping the Python API.

```cython
# solver_fast.pyx
cimport numpy as np
import numpy as np
from cpython cimport bool

cdef class FastSolver:
    cdef:
        double[:] _Etot, _DEtot, _sigma, _statev
        double[:,:] _Lt, _F0, _F1, _DR
        double[:] _Wm, _props
        object _scc  # simcoon._core module

    def __init__(self, int nstatev, props):
        # Pre-allocate all arrays as memoryviews
        self._Etot = np.zeros(6, dtype=np.float64)
        self._DEtot = np.zeros(6, dtype=np.float64)
        self._sigma = np.zeros(6, dtype=np.float64)
        self._statev = np.zeros(nstatev, dtype=np.float64)
        self._Lt = np.zeros((6, 6), dtype=np.float64, order='F')
        self._F0 = np.eye(3, dtype=np.float64, order='F')
        self._F1 = np.eye(3, dtype=np.float64, order='F')
        self._DR = np.eye(3, dtype=np.float64, order='F')
        self._Wm = np.zeros(4, dtype=np.float64)
        self._props = np.asarray(props, dtype=np.float64)

        import simcoon._core as scc
        self._scc = scc

    cpdef void solve_step_fast(self, str umat_name, double[:] DEtot_target,
                                int ninc, double time_step):
        cdef:
            int i, j
            double Dtinc = 1.0 / ninc
            double DTime = Dtinc * time_step
            double Time = 0.0

        for i in range(ninc):
            # Apply strain increment (typed memoryview operations)
            for j in range(6):
                self._DEtot[j] = DEtot_target[j] * Dtinc

            # Reset DR to identity
            for j in range(3):
                self._DR[j, j] = 1.0

            # Call UMAT (still Python call but minimal overhead)
            self._call_umat_fast(umat_name, Time, DTime)

            # Update total strain
            for j in range(6):
                self._Etot[j] += self._DEtot[j]

            Time += DTime

    cdef void _call_umat_fast(self, str umat_name, double Time, double DTime):
        # Reshape to 2D for binding (views, no copy)
        cdef np.ndarray etot_2d = np.asarray(self._Etot).reshape(6, 1, order='F')
        cdef np.ndarray Detot_2d = np.asarray(self._DEtot).reshape(6, 1, order='F')
        # ... etc

        self._scc.umat_inplace(umat_name, etot_2d, Detot_2d, ...)
```

**Advantages:**
- 10-50x faster loop execution
- Maintains Python API compatibility
- Can call existing C++ bindings
- Gradual migration possible

**Disadvantages:**
- Build complexity (requires Cython compilation)
- Debugging harder
- Need to maintain both Python and Cython versions

---

### Strategy B: NumPy Vectorization (Batch Increments)

**Effort: Low | Impact: Medium (1.5-2x speedup)**

Process multiple increments in a single C++ call by pre-computing the full strain path.

```python
class Solver:
    def solve_vectorized(self, sv_init=None):
        """Vectorized solver for strain-controlled loading."""
        # Pre-compute all strain increments
        all_DEtot = []
        all_times = []

        for block in self.blocks:
            for step in block.steps:
                ninc = step.Dn_init
                for i in range(ninc):
                    Dtinc = 1.0 / ninc
                    all_DEtot.append(step.DEtot_end * Dtinc)
                    all_times.append(step.time * Dtinc)

        # Stack into batch arrays
        n_total = len(all_DEtot)
        etot_batch = np.zeros((6, n_total), order='F')
        Detot_batch = np.column_stack(all_DEtot).astype(np.float64, order='F')
        sigma_batch = np.zeros((6, n_total), order='F')
        # ... other arrays

        # Cumulative strain
        etot_batch = np.cumsum(Detot_batch, axis=1)

        # Single batched UMAT call (if UMAT supports it)
        # OR: Loop in C++ instead of Python
        for i in range(n_total):
            scc.umat_inplace(...)  # Still per-point, but arrays pre-allocated

        return self._extract_history(sigma_batch, etot_batch, ...)
```

**Advantages:**
- Pure Python, no compilation
- Reduces Python loop overhead
- Works with existing bindings

**Disadvantages:**
- Only works for fully strain-controlled
- Memory overhead for large simulations
- Doesn't help with Newton-Raphson iterations

---

### Strategy C: Numba JIT Compilation

**Effort: Low-Medium | Impact: Medium-High (2-4x speedup)**

Use Numba to JIT-compile the solver loop.

```python
from numba import jit, float64
from numba.experimental import jitclass

# Define state as a Numba-compatible structure
state_spec = [
    ('Etot', float64[:]),
    ('DEtot', float64[:]),
    ('sigma', float64[:]),
    ('statev', float64[:]),
    ('Lt', float64[:,:]),
    ('Wm', float64[:]),
]

@jitclass(state_spec)
class FastState:
    def __init__(self, nstatev):
        self.Etot = np.zeros(6)
        self.DEtot = np.zeros(6)
        self.sigma = np.zeros(6)
        self.statev = np.zeros(nstatev)
        self.Lt = np.zeros((6, 6))
        self.Wm = np.zeros(4)

@jit(nopython=False)  # object mode to call C++ binding
def solve_step_numba(state, umat_func, props, DEtot_target, ninc, time_step):
    Dtinc = 1.0 / ninc
    DTime = Dtinc * time_step
    Time = 0.0

    history_Etot = np.zeros((ninc + 1, 6))
    history_sigma = np.zeros((ninc + 1, 6))
    history_Etot[0] = state.Etot.copy()
    history_sigma[0] = state.sigma.copy()

    for i in range(ninc):
        # Apply increment
        state.DEtot[:] = DEtot_target * Dtinc

        # Call UMAT (escapes to Python/C++)
        umat_func(state, props, Time, DTime)

        # Update
        state.Etot += state.DEtot
        Time += DTime

        # Store history
        history_Etot[i + 1] = state.Etot.copy()
        history_sigma[i + 1] = state.sigma.copy()

    return history_Etot, history_sigma
```

**Advantages:**
- Easy to implement
- No compilation step (JIT at runtime)
- Can fall back to object mode for C++ calls

**Disadvantages:**
- First call has JIT overhead
- Object mode for C++ calls limits speedup
- Numba updates can break code

---

### Strategy D: C++ Solver with Python Configuration

**Effort: High | Impact: Very High (5-10x speedup)**

Move the entire solver loop to C++, configure from Python.

```cpp
// solver_loop.cpp
void solve_loop_inplace(
    const std::string& umat_name,
    py::array_t<double>& etot_history,      // (6, n_increments) output
    py::array_t<double>& sigma_history,     // (6, n_increments) output
    const py::array_t<double>& DEtot_total, // (6,) total strain
    const py::array_t<double>& props,
    int nstatev,
    int n_increments,
    double total_time
) {
    // All arrays as views
    auto etot_hist = carma::arr_to_mat_view(etot_history);
    auto sigma_hist = carma::arr_to_mat_view(sigma_history);
    auto DEtot = carma::arr_to_col_view(DEtot_total);
    auto props_vec = carma::arr_to_col_view(props);

    // Local state
    vec etot = zeros(6);
    vec Detot = zeros(6);
    vec sigma = zeros(6);
    vec statev = zeros(nstatev);
    mat Lt = zeros(6, 6);
    mat DR = eye(3, 3);
    vec Wm = zeros(4);

    double Dtinc = 1.0 / n_increments;
    double DTime = total_time / n_increments;
    double Time = 0.0;

    // Store initial state
    etot_hist.col(0) = etot;
    sigma_hist.col(0) = sigma;

    for (int i = 0; i < n_increments; i++) {
        // Apply increment
        Detot = DEtot * Dtinc;

        // Call UMAT directly (no binding overhead)
        simcoon::umat_elasticity_iso(etot, Detot, sigma, Lt, L, sigma_in,
                                      DR, nprops, props_vec, nstatev, statev,
                                      0.0, 0.0, Time, DTime,
                                      Wm(0), Wm(1), Wm(2), Wm(3),
                                      3, 3, i == 0, 0, tnew_dt);

        // Update
        etot += Detot;
        Time += DTime;

        // Store history
        etot_hist.col(i + 1) = etot;
        sigma_hist.col(i + 1) = sigma;
    }
}
```

**Python wrapper:**
```python
def solve_fast(umat_name, props, DEtot_total, n_increments, total_time, nstatev=1):
    """Fast solver for strain-controlled loading."""
    # Pre-allocate history arrays
    etot_history = np.zeros((6, n_increments + 1), order='F')
    sigma_history = np.zeros((6, n_increments + 1), order='F')

    # Call C++ solver loop
    scc.solve_loop_inplace(
        umat_name, etot_history, sigma_history,
        DEtot_total, props, nstatev, n_increments, total_time
    )

    return etot_history, sigma_history
```

**Advantages:**
- Maximum performance (close to pure C++)
- No Python overhead in loop
- Still configurable from Python

**Disadvantages:**
- Significant C++ development
- Less flexible (need C++ changes for new features)
- Harder to debug

---

### Strategy E: Reduce History Storage Overhead

**Effort: Low | Impact: Low-Medium (10-20% speedup)**

Store history less frequently or use pre-allocated arrays.

```python
class Solver:
    def __init__(self, ..., history_interval=1):
        self.history_interval = history_interval
        # Pre-allocate history arrays
        self._max_history = 10000
        self._history_Etot = np.zeros((self._max_history, 6))
        self._history_sigma = np.zeros((self._max_history, 6))
        self._history_count = 0

    def _store_history(self, sv):
        """Store to pre-allocated arrays instead of creating objects."""
        if self._history_count < self._max_history:
            np.copyto(self._history_Etot[self._history_count], sv.Etot)
            np.copyto(self._history_sigma[self._history_count], sv.sigma)
            self._history_count += 1

    def solve(self, ...):
        # ...
        for i, increment in enumerate(increments):
            # ... solve increment ...

            # Store every N increments
            if i % self.history_interval == 0:
                self._store_history(sv)
```

**Advantages:**
- Very easy to implement
- No dependencies
- Reduces memory allocations

**Disadvantages:**
- Limited speedup
- May lose intermediate history points

---

### Strategy F: PyPy Compatibility

**Effort: Medium | Impact: High (3-5x speedup)**

Make the solver compatible with PyPy for faster execution.

```python
# Avoid features PyPy doesn't optimize well:
# - Avoid __slots__ in dataclasses
# - Use simple loops instead of numpy for small arrays
# - Minimize object creation in hot path

class FastStateVariables:
    """PyPy-friendly state variables."""
    def __init__(self, nstatev):
        # Use lists instead of numpy for small fixed-size arrays
        self.Etot = [0.0] * 6
        self.DEtot = [0.0] * 6
        self.sigma = [0.0] * 6
        self.statev = [0.0] * nstatev
        # ... etc

    def apply_increment(self, DEtot_target, Dtinc):
        for i in range(6):
            self.DEtot[i] = DEtot_target[i] * Dtinc
```

**Advantages:**
- Significant speedup for pure Python code
- No compilation step
- JIT improves over time

**Disadvantages:**
- NumPy operations may be slower
- C extension compatibility issues
- Need to test with PyPy

---

## 4. Recommendation

### Short-term (Low effort, Quick wins)

1. **Strategy E: Pre-allocated history** - 10-20% speedup
2. **Strategy B: Batch processing** for strain-controlled cases - 1.5x speedup

### Medium-term (Recommended)

3. **Strategy C: Numba JIT** - 2-4x speedup with minimal code changes
4. **Strategy A: Cython hot path** - Best balance of performance and flexibility

### Long-term (Maximum performance)

5. **Strategy D: C++ solver loop** - For production/HPC use cases

---

## 5. Implementation Priority

```
Phase 1 (Now):
├── Pre-allocated history arrays
├── Reduce attribute access in hot path
└── Profile to identify remaining bottlenecks

Phase 2 (Next):
├── Numba JIT for solver loop
├── Cython for _solve_increment and _call_umat
└── Benchmark against C++ solver

Phase 3 (Future):
├── C++ solver loop option
├── Multi-threaded batch processing
└── GPU acceleration for large batches
```

---

## 6. Quick Wins (Implement Now)

### 5.1 Reduce Attribute Access

```python
# Before (slow - attribute access in loop)
def _solve_increment(self, ...):
    sv.DEtot[:] = DEtot_target * Dtinc
    sv.DT = Dtinc * DT_target

# After (faster - local variables)
def _solve_increment(self, ...):
    DEtot = sv.DEtot  # Cache reference
    DEtot[:] = DEtot_target
    DEtot *= Dtinc
```

### 5.2 Avoid Repeated Method Lookups

```python
# Before
for i in range(n):
    np.copyto(dst, src)

# After
_copyto = np.copyto  # Cache function reference
for i in range(n):
    _copyto(dst, src)
```

### 5.3 Use __slots__ in Hot Path Classes

```python
@dataclass
class HistoryPoint:
    __slots__ = ['Etot', 'sigma', 'Wm', 'statev', 'R', 'T']
    Etot: np.ndarray
    sigma: np.ndarray
    # ...
```

---

## 7. Benchmarking Script

```python
"""Benchmark different solver implementations."""
import numpy as np
import time

def benchmark_solver(solver_class, name, n_runs=10):
    props = np.array([210000.0, 0.3, 0.0])
    step = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
                    Dn_init=1000, Dn_mini=1000, Dn_inc=1000)
    block = Block(steps=[step], umat_name='ELISO', props=props, nstatev=1)

    times = []
    for _ in range(n_runs):
        solver = solver_class(blocks=[block])
        start = time.perf_counter()
        history = solver.solve()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    avg = np.mean(times) * 1000
    std = np.std(times) * 1000
    per_inc = avg / 1000 * 1000  # µs per increment

    print(f'{name:20s}: {avg:6.2f} ± {std:4.2f} ms ({per_inc:.1f} µs/inc)')

# Run benchmarks
benchmark_solver(Solver, 'Current Python')
# benchmark_solver(CythonSolver, 'Cython')
# benchmark_solver(NumbaSolver, 'Numba')
```

---

## 8. Conclusion

### Key Findings

**Current Performance:**
- Python solver: 14.5 µs/increment
- C++ solver: ~12 µs/increment
- Ratio: 1.2x (Python is 20% slower than C++)

**Bottleneck Distribution:**
| Source | Time | Actionable? |
|--------|------|-------------|
| C++ UMAT + binding | 4.8 µs | No (already optimized) |
| np.linalg.solve | 3.3 µs | No (NumPy overhead) |
| Python loop | 4.6 µs | **Yes** (main target) |
| History storage | 1.3 µs | Yes (minor) |
| Array reshaping | 0.5 µs | No (negligible) |

**Migration Assessment:**
- **Eigen migration: NOT RECOMMENDED** - 1900-2950 hours effort for <0.1 µs gain
- **Armadillo+Carma: Keep** - Already provides zero-copy data transfer
- **Nanobind: Future option** - Could reduce binding overhead by 0.1-0.2 µs

### Optimization Path

The Python solver can be optimized from 14.5 µs to 10-11 µs per increment through:

1. **Quick wins** (Section 6): Pre-allocated arrays, cached references, `__slots__`
2. **Numba/Cython** (Strategies A/C): JIT compilation of hot path
3. **C++ loop** (Strategy D): Move entire loop to C++ for maximum performance

The recommended approach is **Cython for the hot path** (Strategy A), which provides:
- 2-3x speedup on Python loop portion
- Maintains Python API flexibility
- Gradual migration possible
- Good debugging support

**Expected outcome:** ~10-11 µs/increment, achieving **1.0-1.1x ratio vs C++**

### Architectural Recommendations

1. **Keep current Armadillo+Carma+pybind11 stack** - It's well-optimized
2. **Focus optimization effort on Python loop**, not C++ binding layer
3. **Cython is the best balance** of effort vs performance gain
4. **C++ solver loop** only if HPC requirements demand it
5. **Monitor nanobind development** as potential future pybind11 replacement

---

## Appendix A: Codebase Statistics

| Metric | Value |
|--------|-------|
| Total C++ files | 1,345 |
| Files using Armadillo | 226 (16.8%) |
| `arma::mat` occurrences | 1,267 |
| `arma::vec` occurrences | 919 |
| UMAT models | 41 |
| Continuum mechanics files | 160+ |

---

## Appendix B: Profiling Methodology

Measurements taken using:
- `timeit` module with 10,000+ iterations
- Isolated micro-operations to measure individual components
- CPython 3.13 on macOS Darwin 23.0.0
- Apple Silicon (arm64)

```python
# Example profiling code
import timeit

setup = '''
import numpy as np
import simcoon._core as scc
# ... setup arrays ...
'''

# Measure umat_inplace
time_umat = timeit.timeit(
    'scc.umat_inplace(...)',
    setup=setup,
    number=10000
) / 10000 * 1e6  # µs
```
