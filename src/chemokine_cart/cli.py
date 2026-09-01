"""Command-line interface for reproducible chemokine–CAR-T analyses."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from .io import (
    DataValidationError,
    load_config,
    make_demo_dataset,
    resolve_project_path,
    validate_configured_inputs,
)


def _add_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="YAML configuration path (default: config/config.yaml)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chemokine-cart",
        description="Validate and reproduce chemokine–CAR-T analyses.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate", help="validate config, manifest, and canonical inputs"
    )
    _add_config_argument(validate_parser)
    validate_parser.add_argument(
        "--manifest",
        default=None,
        help="project-relative manifest override",
    )
    validate_parser.add_argument(
        "--check-files",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="check local files and their SHA-256 digests",
    )

    demo_parser = subparsers.add_parser(
        "make-demo",
        help="create deterministic synthetic software-test data",
    )
    _add_config_argument(demo_parser)
    demo_parser.add_argument("--output-dir", default="data/demo", help="project-relative output")
    demo_parser.add_argument("--seed", type=int, default=None, help="deterministic random seed")
    demo_parser.add_argument("--force", action="store_true", help="replace existing demo files")

    run_parser = subparsers.add_parser("run", help="validate inputs and execute the pipeline")
    _add_config_argument(run_parser)
    run_parser.add_argument("--manifest", default=None, help="project-relative manifest override")
    run_parser.add_argument("--output-dir", default=None, help="project-relative output override")
    run_parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="override project.random_seed for this run",
    )
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and print planned work without producing results",
    )
    return parser


def _config_path(raw: str) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = Path.cwd() / path
    return path.resolve()


def _cmd_validate(args: argparse.Namespace) -> int:
    config = load_config(_config_path(args.config))
    report = validate_configured_inputs(
        config,
        manifest_override=args.manifest,
        check_files=args.check_files,
    )
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0


def _cmd_make_demo(args: argparse.Namespace) -> int:
    config = load_config(_config_path(args.config))
    seed = args.seed
    if seed is None:
        seed = int(config.get("project", {}).get("random_seed", 20260901))
    manifest = make_demo_dataset(
        Path(config["_project_root"]),
        args.output_dir,
        seed=seed,
        force=args.force,
    )
    relative = manifest.relative_to(Path(config["_project_root"]))
    print(f"Synthetic software-test manifest: {relative.as_posix()}")
    print("Evidence eligibility: false")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    config_path = _config_path(args.config)
    config = load_config(config_path)
    validate_configured_inputs(
        config,
        manifest_override=args.manifest,
        check_files=True,
    )
    root = Path(config["_project_root"])
    manifest_value = args.manifest or str(config["manifest"]["path"])
    output_value = args.output_dir or str(config.get("output", {}).get("root", "results"))
    manifest_path = resolve_project_path(root, manifest_value, must_exist=True)
    output_dir = resolve_project_path(root, output_value)

    from .pipeline import run_pipeline

    run_pipeline(
        config_path=config_path,
        manifest_path=manifest_path,
        output_dir=output_dir,
        dry_run=args.dry_run,
        random_seed=args.random_seed,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _cmd_validate(args)
        if args.command == "make-demo":
            return _cmd_make_demo(args)
        if args.command == "run":
            return _cmd_run(args)
    except (
        DataValidationError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
    ) as exc:
        parser.exit(2, f"error: {exc}\n")
    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
