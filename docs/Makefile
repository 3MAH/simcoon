# Makefile for Simcoon documentation (Sphinx + Doxygen/Breathe)
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Doxygen configuration
DOXYFILE      = Doxyfile
DOXYGEN       ?= doxygen
DOXYGENDIR    = $(BUILDDIR)/doxygen
DOXYGENXML    = $(DOXYGENDIR)/xml

# Put it first so that "make" without argument is like "make help".
help:
	@echo "Simcoon Documentation Build System"
	@echo "==================================="
	@echo ""
	@echo "Available targets:"
	@echo "  html        - Build HTML documentation (runs doxygen first if needed)"
	@echo "  doxygen     - Generate Doxygen XML for C++ API"
	@echo "  clean       - Remove all build artifacts"
	@echo "  clean-doxy  - Remove only Doxygen output"
	@echo "  rebuild     - Clean and rebuild everything"
	@echo "  serve       - Start local server to view docs (port 11000)"
	@echo "  open        - Open documentation in browser (macOS)"
	@echo ""
	@echo "Advanced targets:"
	@echo "  html-only   - Build HTML without regenerating Doxygen"
	@echo "  linkcheck   - Check for broken links"
	@echo "  doctest     - Run doctests in documentation"
	@echo ""
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile doxygen clean-doxy rebuild serve open html-only clean html linkcheck doctest serve-html

# Generate Doxygen XML documentation for Breathe
doxygen:
	@echo "Running Doxygen to generate XML for Breathe..."
	@mkdir -p $(DOXYGENDIR)
	$(DOXYGEN) $(DOXYFILE)
	@echo "Doxygen XML generated in $(DOXYGENXML)"

# Build HTML with automatic Doxygen generation
html: doxygen
	@echo "Building Sphinx HTML documentation..."
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo ""
	@echo "Documentation built successfully!"
	@echo "Open $(BUILDDIR)/html/index.html to view."

# Build HTML without regenerating Doxygen (faster for content-only changes)
html-only:
	@echo "Building Sphinx HTML (skipping Doxygen)..."
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Clean only Doxygen output
clean-doxy:
	@echo "Removing Doxygen output..."
	rm -rf $(DOXYGENDIR)

# Full clean
clean:
	@echo "Cleaning all build artifacts..."
	rm -rf $(BUILDDIR)

# Rebuild everything from scratch
rebuild: clean html
	@echo "Full rebuild complete!"

# Spin up a local http server to view the rendered documentation.
# This is required for interactive examples to work.
serve: html
	@echo "Starting local server at http://localhost:11000"
	@echo "Press Ctrl+C to stop."
	python -m http.server 11000 --directory "$(BUILDDIR)"/html

# Alias for backward compatibility
serve-html: serve

# Open documentation in browser (macOS)
open: html
	@open "$(BUILDDIR)/html/index.html"

# Check for broken links
linkcheck: doxygen
	@$(SPHINXBUILD) -M linkcheck "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Run doctests
doctest: doxygen
	@$(SPHINXBUILD) -M doctest "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
# Note: doxygen dependency ensures C++ docs are generated first
%: Makefile
	@$(MAKE) doxygen --no-print-directory
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
