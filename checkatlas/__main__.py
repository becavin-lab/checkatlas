import argparse  # pragma: no cover

from . import atlas # pragma: no cover
from . import multiqc # pragma: no cover


def main() -> None:  # pragma: no cover
    """
    The main function executes on commands:
    `python -m checkatlas` and `$ checkatlas `.

    This is your program's entry point.

    You can change this function to do whatever you want.
    Examples:
        * Run a test suite
        * Run a server
        * Do some other stuff
        * Run a command line application (Click, Typer, ArgParse)
        * List all available tasks
        * Run an application (Flask, FastAPI, Django, etc.)
    """
    parser = argparse.ArgumentParser(
        description="checkatlas.",
        epilog="Enjoy the checkatlas functionality!",
    )
    # This is required positional argument
    parser.add_argument(
        "path",
        type=str,
        help="Path containing Scanpy and Seurat atlases",
        default="becavin-lab",
    )
    # This is optional named argument
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Optionally adds verbosity",
    )
    args = parser.parse_args()
    print(f"Checking your single-cell atlases in {args.path}!")
    if args.verbose:
        print("Verbose mode is on.")

    print("Searching Seurat and Scanpy files")
    atlas_list = atlas.list_atlases(args.path)
    print('Found', len(atlas_list),'atlas(es)',atlas_list)

    print("Convert all Seurat atlases")
    atlas.convert_seurat_atlases(args.path, atlas_list)

    print("Clean Scanpy atlases")
    #atlas.clean_scanpy(args.path, atlas_list)

    print("Parse information in Atlas")
    atlas.parse_atlases(args.path, atlas_list)

    print("Run metric calculation")


    print("Run MultiQC")
    #multiqc.run_multiqc(args.path)


if __name__ == "__main__":  # pragma: no cover
    main()
