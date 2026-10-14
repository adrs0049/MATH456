# Local build. Deployment is the GitHub workflow in .github/workflows.
#
# The site uses the myst-theme fork in the myst-theme submodule (as MATH725).
# Build it once after cloning:  git submodule update --init && make theme
# That assembles a production copy in _theme/ whose server is `node server.js`
# (the submodule's own template.yml starts `remix dev`, which rebuilds the
# assets while myst copies them and breaks local builds).
#
# Every target executes code cells (myst's --execute) against a throwaway
# Jupyter server, because the figures are code cells and
# mystmd has no config switch for execution, only the CLI flag. Executed
# outputs are cached under _build/execute and reused until a cell changes.
# The *-noexec targets skip execution and leave those figures out.
#
# Serves on http://localhost:3002 by default, so it can run beside the
# private MATH456 (3000) and MATH690CB (3001) books. Override with:
#   make start PORT=8080 SERVER_PORT=8180
PORT        ?= 3002
SERVER_PORT ?= 3102
# Port for the throwaway Jupyter server that executes notebooks.
JPORT       ?= 8913
JTOKEN      ?= mysttoken
JLOG        := /tmp/myst-jupyter-$(JPORT).log
# Pages executed at once. mystmd defaults to one per core, which on a big
# machine starves the slow figure cells (FitzHugh-Nagumo, saddle-node) until
# they time out.
EXEC_PARALLEL ?= 8

.PHONY: theme start build pdf pdf-page start-noexec build-noexec pdf-noexec pdf-page-noexec \
        build-exec pdf-exec pdf-page-exec jupyter-start jupyter-stop clean

start:            ## live preview at http://localhost:$(PORT), code cells executed
	@$(MAKE) --no-print-directory HOST=127.0.0.1 _exec MYST_CMD="start --port $(PORT) --server-port $(SERVER_PORT) --execute --execute-parallel $(EXEC_PARALLEL)"

build:            ## static HTML into _build/, code cells executed
	@$(MAKE) --no-print-directory _exec MYST_CMD="build --html --execute --execute-parallel $(EXEC_PARALLEL)"

pdf:              ## every declared per-page PDF export, code cells executed
	@$(MAKE) --no-print-directory _exec MYST_CMD="build --pdf --execute --execute-parallel $(EXEC_PARALLEL)"

pdf-page:         ## one page only:  make pdf-page FILE=overdamped/singular-limit.md
	@$(MAKE) --no-print-directory _exec MYST_CMD="build $(FILE) --pdf --execute --execute-parallel $(EXEC_PARALLEL)"

# Old names, kept so existing habits still work.
build-exec: build
pdf-exec: pdf
pdf-page-exec: pdf-page

start-noexec:     ## live preview without executing code cells
	myst start --port $(PORT) --server-port $(SERVER_PORT)

build-noexec:     ## static HTML without executing code cells
	myst build --html

pdf-noexec:       ## PDF exports without executing code cells
	myst build --pdf

pdf-page-noexec:  ## one page, no execution:  make pdf-page-noexec FILE=...
	myst build $(FILE) --pdf

# Run one myst command with the Jupyter server up, and take the server down
# afterwards whatever happens (including Ctrl-C out of `make start`).
# Do NOT use `pkill -f` on the jupyter pattern: it matches the invoking shell.
_exec:
	@$(MAKE) --no-print-directory jupyter-start
	@trap '$(MAKE) --no-print-directory jupyter-stop' EXIT INT TERM; \
	  JUPYTER_BASE_URL=http://127.0.0.1:$(JPORT) JUPYTER_TOKEN=$(JTOKEN) \
	  myst $(MYST_CMD)

jupyter-start:    ## start the throwaway Jupyter server on :$(JPORT)
	@if curl -sf -o /dev/null http://127.0.0.1:$(JPORT)/api; then \
	  echo ">> Jupyter server already up on :$(JPORT)"; \
	else \
	  echo ">> starting Jupyter server on :$(JPORT) for notebook execution"; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
	  jupyter server --no-browser --port $(JPORT) \
	    --IdentityProvider.token=$(JTOKEN) --ServerApp.disable_check_xsrf=True \
	    >$(JLOG) 2>&1 & \
	  for i in $$(seq 1 60); do \
	    curl -sf -o /dev/null http://127.0.0.1:$(JPORT)/api && break; sleep 0.5; \
	  done; \
	fi

jupyter-stop:     ## stop the throwaway Jupyter server
	@jupyter server stop $(JPORT) >/dev/null 2>&1 || true
	@rm -f /tmp/myst-jupyter-$(JPORT).pid

theme:            ## build the myst-theme submodule into _theme/ (once, and after updating it)
	cd myst-theme && npm ci && npm run build
	cd myst-theme/themes/book && rm -rf public build && npm run prod:build
	rm -rf _theme && mkdir -p _theme
	cp -r myst-theme/themes/book/public myst-theme/themes/book/build \
	      myst-theme/themes/book/template.yml myst-theme/template/package.json \
	      myst-theme/template/server.js _theme/
	sed -i 's/"@myst-theme\/template"/"@myst-theme\/book"/; s/"VERSION"/"0.0.0"/' _theme/package.json
	cd _theme && npm install --ignore-scripts

clean:
	myst clean --all --yes
