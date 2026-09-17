from __future__ import annotations

from itertools import product
import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _as_arg_names(raw_names):
    if isinstance(raw_names, str):
        return [name.strip() for name in raw_names.split(",") if name.strip()]
    return list(raw_names)


def _as_value_tuple(raw_value, n_args):
    values = getattr(raw_value, "values", raw_value)
    if n_args == 1:
        return (values,)
    return tuple(values)


def _parametrize_cases(func):
    marks = [
        mark
        for mark in getattr(func, "pytestmark", [])
        if getattr(mark, "name", None) == "parametrize"
    ]
    if not marks:
        return [(func.__name__, (), {})]

    case_groups = []
    for mark in marks:
        arg_names = _as_arg_names(mark.args[0])
        arg_values = mark.args[1]
        group = []
        for index, raw_value in enumerate(arg_values):
            value_tuple = _as_value_tuple(raw_value, len(arg_names))
            case_id = getattr(raw_value, "id", None) or str(index)
            group.append((case_id, dict(zip(arg_names, value_tuple))))
        case_groups.append(group)

    cases = []
    for combo in product(*case_groups):
        label = func.__name__ + "[" + "-".join(item[0] for item in combo) + "]"
        kwargs = {}
        for _, group_kwargs in combo:
            kwargs.update(group_kwargs)
        cases.append((label, (), kwargs))
    return cases


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
                test_func = getattr(module, name)
                for case_label, args, kwargs in _parametrize_cases(test_func):
                    try:
                        test_func(*args, **kwargs)
                        print(f"PASS {path.name}:{case_label}")
                    except Exception as exc:
                        failures.append((path.name, case_label, exc))
                        print(f"FAIL {path.name}:{case_label}: {type(exc).__name__}: {exc}")
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
