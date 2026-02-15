# ribctl/lib/npet2/__main__.py
"""
npet2 CLI entry point.

Usage:
    python -m ribctl.lib.npet2 run 7K00 4UG0 --workers 4
    python -m ribctl.lib.npet2 run --from-file structures.txt --output-dir ./results
    python -m ribctl.lib.npet2 run 7K00 --cylinder-radius 40 --voxel-size 0.5
    python -m ribctl.lib.npet2 run 7K00 --mmcif /path/to/7K00.cif --profile /path/to/profile.json --landmarks /path/to/landmarks.json
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Optional


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="npet2",
        description="Ribosome exit tunnel geometry pipeline",
    )
    sub = p.add_subparsers(dest="command")

    # --- run ---
    run_p = sub.add_parser("run", help="Run the pipeline on one or more structures")

    # Input selection
    run_p.add_argument("rcsb_ids", nargs="*", help="RCSB IDs to process")
    run_p.add_argument("--from-file", type=str, default=None,
                       help="Read RCSB IDs from a text file (one per line)")

    # Provider mode
    run_p.add_argument("--mode", choices=["riboxyz", "standalone"], default="riboxyz",
                       help="riboxyz: use local riboxyz assets. standalone: use mmcif + profile files or API.")
    run_p.add_argument("--mmcif", type=str, default=None,
                       help="(standalone) Path to mmCIF file. For multiple structures, use a directory.")
    run_p.add_argument("--profile", type=str, default=None,
                       help="(standalone) Path to profile JSON file or directory of profiles.")
    run_p.add_argument("--landmarks", type=str, default=None,
                       help="(standalone) Path to landmarks JSON file or directory.")
    run_p.add_argument("--api-url", type=str, default=None,
                       help="(standalone) riboxyz API base URL for fetching profiles/landmarks.")

    # Output
    run_p.add_argument("--output-dir", type=str, default=None,
                       help="Root output directory (default: NPET2_RUNS_ROOT)")

    # Parallelism
    run_p.add_argument("--workers", "-j", type=int, default=1,
                       help="Number of parallel workers (default: 1)")

    # Config overrides (the commonly-tuned ones)
    cfg = run_p.add_argument_group("config overrides")
    cfg.add_argument("--cylinder-radius", type=float, default=None)
    cfg.add_argument("--cylinder-height", type=float, default=None)
    cfg.add_argument("--ptc-extension", type=float, default=None)
    cfg.add_argument("--voxel-size", type=float, default=None,
                     help="Level-0 grid voxel size in Angstroms")
    cfg.add_argument("--refine-voxel-size", type=float, default=None)
    cfg.add_argument("--no-mesh", action="store_true", help="Disable mesh generation")
    cfg.add_argument("--no-refine", action="store_true", help="Skip Stage55 grid refinement")
    cfg.add_argument("--dbscan-coarse-eps", type=float, default=None)
    cfg.add_argument("--dbscan-coarse-min-samples", type=int, default=None)
    cfg.add_argument("--dbscan-refine-eps", type=float, default=None)
    cfg.add_argument("--dbscan-refine-min-samples", type=int, default=None)
    cfg.add_argument("--config-json", type=str, default=None,
                     help="Path to a full RunConfig JSON (overrides individual flags)")

    # --- show-config ---
    sub.add_parser("show-config", help="Print the default RunConfig as JSON")

    return p


def _build_config(args) -> "RunConfig":
    from ribctl.lib.npet2.core.config import RunConfig, GridLevelConfig

    if args.config_json:
        import json
        data = json.loads(Path(args.config_json).read_text())
        # Reconstruct GridLevelConfig objects
        if "grid_levels" in data:
            data["grid_levels"] = [GridLevelConfig(**gl) for gl in data["grid_levels"]]
        return RunConfig(**data)

    kwargs = {}
    if args.cylinder_radius is not None:
        kwargs["cylinder_radius_A"] = args.cylinder_radius
    if args.cylinder_height is not None:
        kwargs["cylinder_height_A"] = args.cylinder_height
    if args.ptc_extension is not None:
        kwargs["cylinder_ptc_extension_A"] = args.ptc_extension
    if args.refine_voxel_size is not None:
        kwargs["refine_voxel_size_A"] = args.refine_voxel_size
    if args.no_mesh:
        kwargs["mesh_level0_enable"] = False
        kwargs["mesh_level1_enable"] = False
    if args.dbscan_coarse_eps is not None:
        kwargs["dbscan_level0_coarse_eps_A"] = args.dbscan_coarse_eps
    if args.dbscan_coarse_min_samples is not None:
        kwargs["dbscan_level0_coarse_min_samples"] = args.dbscan_coarse_min_samples
    if args.dbscan_refine_eps is not None:
        kwargs["dbscan_level0_refine_eps_A"] = args.dbscan_refine_eps
    if args.dbscan_refine_min_samples is not None:
        kwargs["dbscan_level0_refine_min_samples"] = args.dbscan_refine_min_samples

    if args.voxel_size is not None:
        kwargs["grid_levels"] = [
            GridLevelConfig(name="level_0", voxel_size_A=args.voxel_size,
                            occupancy_backend="legacy_kdtree"),
        ]

    return RunConfig(**kwargs)


def _collect_rcsb_ids(args) -> list[str]:
    ids = list(args.rcsb_ids) if args.rcsb_ids else []
    if args.from_file:
        p = Path(args.from_file)
        if not p.exists():
            print(f"Error: --from-file {p} does not exist", file=sys.stderr)
            sys.exit(1)
        for line in p.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line)
    if not ids:
        print("Error: no RCSB IDs specified", file=sys.stderr)
        sys.exit(1)
    return [x.upper() for x in ids]


def _make_providers(args, rcsb_id: str):
    """Return (structure_provider, landmark_provider) for a given structure."""
    if args.mode == "riboxyz":
        from ribctl.lib.npet2.adapters.riboxyz_providers import (
            RiboxyzStructureProvider,
            RiboxyzLandmarkProvider,
        )
        return RiboxyzStructureProvider(), RiboxyzLandmarkProvider()

    # standalone mode
    from ribctl.lib.npet2.adapters.standalone_providers import (
        FileStructureProvider,
        FileLandmarkProvider,
    )

    api_base = args.api_url

    # Resolve mmcif path
    mmcif_path = None
    if args.mmcif:
        p = Path(args.mmcif)
        if p.is_dir():
            # Look for {RCSB_ID}.cif or {rcsb_id}.cif
            for candidate in [f"{rcsb_id}.cif", f"{rcsb_id.lower()}.cif"]:
                if (p / candidate).exists():
                    mmcif_path = p / candidate
                    break
            if mmcif_path is None:
                raise FileNotFoundError(
                    f"No mmCIF file found for {rcsb_id} in {p}. "
                    f"Expected {rcsb_id}.cif"
                )
        else:
            mmcif_path = p

    if mmcif_path is None:
        raise ValueError(f"--mmcif is required in standalone mode (for {rcsb_id})")

    # Resolve profile path
    profile_path = None
    if args.profile:
        p = Path(args.profile)
        if p.is_dir():
            for candidate in [f"{rcsb_id}_profile.json", f"{rcsb_id}.json"]:
                if (p / candidate).exists():
                    profile_path = p / candidate
                    break
        else:
            profile_path = p

    # Resolve landmarks path
    landmarks_path = None
    if args.landmarks:
        p = Path(args.landmarks)
        if p.is_dir():
            for candidate in [f"{rcsb_id}_landmarks.json", f"{rcsb_id}.json"]:
                if (p / candidate).exists():
                    landmarks_path = p / candidate
                    break
        else:
            landmarks_path = p

    return (
        FileStructureProvider(mmcif_path, profile_path=profile_path, api_base=api_base),
        FileLandmarkProvider(landmarks_path=landmarks_path, api_base=api_base),
    )


def _run_single(
    rcsb_id: str,
    args,
    config: "RunConfig",
    output_root: Optional[Path],
) -> dict:
    """Run pipeline for a single structure. Returns a result dict."""
    from ribctl.lib.npet2.run import run_npet2

    try:
        sp, lp = _make_providers(args, rcsb_id)

        # Allow output dir override
        if output_root:
            import ribctl.lib.npet2.core.settings as settings
            settings.NPET2_RUNS_ROOT = output_root

        ctx = run_npet2(rcsb_id, config, structure_provider=sp, landmark_provider=lp)
        return {
            "rcsb_id": rcsb_id,
            "status": "success",
            "run_dir": str(ctx.store.run_dir),
        }
    except Exception as e:
        return {
            "rcsb_id": rcsb_id,
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


def _run_worker(packed_args: tuple) -> dict:
    """Wrapper for ProcessPoolExecutor."""
    rcsb_id, args_ns, config_dict, output_root_str = packed_args
    from ribctl.lib.npet2.core.config import RunConfig, GridLevelConfig

    # Reconstruct config from dict
    if "grid_levels" in config_dict:
        config_dict["grid_levels"] = [GridLevelConfig(**gl) for gl in config_dict["grid_levels"]]
    config = RunConfig(**config_dict)

    output_root = Path(output_root_str) if output_root_str else None
    return _run_single(rcsb_id, args_ns, config, output_root)


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "show-config":
        from ribctl.lib.npet2.core.config import RunConfig
        cfg = RunConfig()
        print(json.dumps(asdict(cfg), indent=2))
        return

    if args.command == "run":
        rcsb_ids = _collect_rcsb_ids(args)
        config = _build_config(args)
        output_root = Path(args.output_dir) if args.output_dir else None

        n_workers = min(args.workers, len(rcsb_ids))

        if args.no_refine:
            # We need to modify the pipeline stages -- simplest way is a flag
            # that run.py checks. For now, store it in config as a workaround.
            pass

        print(f"npet2: processing {len(rcsb_ids)} structure(s), workers={n_workers}")

        if n_workers <= 1:
            results = []
            for rid in rcsb_ids:
                r = _run_single(rid, args, config, output_root)
                results.append(r)
                status = r["status"]
                print(f"  {rid}: {status}" + (f" -> {r.get('run_dir', '')}" if status == "success" else f" ({r.get('error', '')})"))
        else:
            # Serialize config for multiprocessing
            config_dict = asdict(config)
            output_root_str = str(output_root) if output_root else None

            packed = [
                (rid, args, config_dict, output_root_str)
                for rid in rcsb_ids
            ]

            results = []
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_worker, p): p[0] for p in packed}
                for fut in as_completed(futures):
                    rid = futures[fut]
                    try:
                        r = fut.result()
                    except Exception as e:
                        r = {"rcsb_id": rid, "status": "failed", "error": str(e)}
                    results.append(r)
                    status = r["status"]
                    print(f"  {rid}: {status}" + (
                        f" -> {r.get('run_dir', '')}" if status == "success"
                        else f" ({r.get('error', '')})"
                    ))

        # Summary
        ok = sum(1 for r in results if r["status"] == "success")
        fail = len(results) - ok
        print(f"\nnpet2: {ok} succeeded, {fail} failed out of {len(results)}")

        if fail > 0:
            sys.exit(1)


if __name__ == "__main__":
    main()