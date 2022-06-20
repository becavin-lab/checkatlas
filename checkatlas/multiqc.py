import os


def run_multiqc(path):
    title = '"Check Single-Cell Atlas"'
    name = "CheckAtlas"
    multiqc_script = (
        "cd " + path + ";"
        "multiqc -i " + title + " -n " + name + " -o CheckAtlas_MultiQC"
        ' -f -e "snippy" --cl-config "ignore_images: false" ' + path
    )
    print(multiqc_script)
    os.system(multiqc_script)
