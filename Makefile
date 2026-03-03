.PHONY: lint format check test clean help

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

# ── Quality ───────────────────────────────────────────────────
lint: ## Run ruff linter
	ruff check bin/ tests/

format: ## Auto-format with ruff
	ruff format bin/ tests/
	ruff check --fix bin/ tests/

check: ## Lint + format check (no changes)
	ruff format --check bin/ tests/
	ruff check bin/ tests/

typecheck: ## Run mypy on bin/
	mypy bin/

# ── Tests ─────────────────────────────────────────────────────
test: ## Run pytest
	pytest

test-cov: ## Run pytest with coverage
	pytest --cov=bin --cov-report=term-missing

# ── Misc ──────────────────────────────────────────────────────
clean: ## Remove caches and build artifacts
	rm -rf .ruff_cache .mypy_cache .pytest_cache __pycache__
	find . -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true
