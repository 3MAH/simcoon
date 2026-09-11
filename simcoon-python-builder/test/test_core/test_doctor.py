"""simcoon.doctor: OpenMP runtime report (structure only; the verdict depends on the env)."""

from simcoon import doctor


def test_loaded_libraries_lists_this_process():
    libs = doctor.loaded_libraries()
    assert isinstance(libs, list)
    if libs:                                   # unsupported platforms return []
        assert any("numpy" in p or "python" in p.lower() or "libSystem" in p for p in libs)


def test_collect_reports_imports_and_openmp_entries():
    report = doctor.collect(modules=("numpy", "simcoon"))
    assert "error" not in report, report.get("error")
    assert [e["module"] for e in report["imports"]] == ["numpy", "simcoon"]
    assert all(e["status"] == "ok" for e in report["imports"])
    assert isinstance(report["openmp"], list)
    for entry in report["openmp"]:
        assert set(entry) == {"path", "brought_by"}
    text = doctor.format_report(report)
    assert "OpenMP runtimes loaded:" in text
    assert ("OK:" in text) or ("PROBLEM:" in text)


def test_format_report_advice_for_duplicates():
    report = {
        "python": "py", "platform": "darwin", "libraries": [],
        "imports": [{"module": "simcoon", "status": "ok", "new_libraries": []},
                    {"module": "torch", "status": "ok", "new_libraries": []}],
        "openmp": [{"path": "/opt/x/libomp.dylib", "brought_by": "simcoon"},
                   {"path": "/env/site-packages/torch/lib/libomp.dylib", "brought_by": "torch"}],
        "blas": [],
    }
    text = doctor.format_report(report)
    assert "PROBLEM" in text and "conda-forge pytorch" in text and "KMP_DUPLICATE_LIB_OK" in text
