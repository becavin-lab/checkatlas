import logging
import os

logger = logging.getLogger("checkatlas")


def run_multiqc(args):
    """
    Run Multiqc software after running checkatlas
    It creates the final html file with al checkatlas
    tables and figures.

    Args:
        args: List of args for checkatlas program
    Returns:
        None: None
    """
    title = '"Check Single-Cell Atlas"'
    name = "CheckAtlas"
    multiqc_script = (
        "cd " + args.path + ";"
        "multiqc -i "
        + title
        + " -n "
        + name
        + " -o "
        + args.multiqc
        + ' -f -e "snippy" --cl-config "ignore_images: false" '
        + args.path
    )
    logger.debug("Run: " + multiqc_script)
    os.system(multiqc_script)
