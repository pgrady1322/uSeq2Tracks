.PHONY: lint format check test clean help

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

# ── Quality ───────────────────────────────────────────────────
lint: ## Run ruff linter
	ruff check src/ tests/

format: ## Auto-format with ruff
	ruff format src/ tests/
	ruff check --fix src/ tests/

check: ## Lint + format check (no changes)
	ruff format --check src/ tests/
	ruff check src/ tests/

typecheck: ## Run mypy on src/
	mypy src/

# ── Tests ─────────────────────────────────────────────────────
test: ## Run pytest
	pytest

test-cov: ## Run pytest with coverage
	pytest --cov=useq2tracks --cov-report=term-missing

# ── Misc ──────────────────────────────────────────────────────
clean: ## Remove caches and build artifacts
	rm -rf .ruff_cache .mypy_cache .pytest_cache __pycache__ dist/ build/ *.egg-info
	find . -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true
