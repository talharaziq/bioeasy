"""Allow BioEasy to run with ``python -m bioeasy``."""

from bioeasy.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
