from __future__ import annotations

import importlib.util
from pathlib import Path


def main() -> int:
    failures = []
    for path in sorted(Path(__file__).parent.glob("test_*.py")):
        spec = importlib.util.spec_from_file_location(path.stem, path)
        module = importlib.util.module_from_spec(spec)
        try:
            assert spec.loader is not None
            spec.loader.exec_module(module)
            for name in sorted(dir(module)):
                if not name.startswith("test_"):
                    continue
                try:
                    getattr(module, name)()
                    print(f"PASS {path.name}:{name}")
                except Exception as exc:
                    failures.append((path.name, name, exc))
                    print(f"FAIL {path.name}:{name}: {type(exc).__name__}: {exc}")
        except Exception as exc:
            failures.append((path.name, "<module>", exc))
            print(f"FAIL {path.name}:<module>: {type(exc).__name__}: {exc}")

    if failures:
        print(f"\n{len(failures)} failure(s)")
        return 1
    print("\nAll tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
