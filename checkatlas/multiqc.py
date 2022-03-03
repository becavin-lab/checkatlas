import os


def run_multiqc(path):
    title = '"Check Single-Cell Atlas"'
    name = "CheckAtlas"
    multiqc_script = (
        "multiqc -i " + title + " -n " + name + " -o MultiQC " + path
    )
    print(multiqc_script)
    os.system(multiqc_script)
