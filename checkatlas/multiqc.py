import os


def run_multiqc(args):
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
    print(multiqc_script)
    os.system(multiqc_script)
