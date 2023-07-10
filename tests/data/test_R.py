import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector


def test_r():
    utils = rpackages.importr("utils")
    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
    # R package names
    packnames = ("Seurat", "SeuratObject")
    # Selectively install what needs to be install.
    # We are fancy, just because we can.
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        rcode = """Sys.getenv("R_LIBS_USER")"""
        print(robjects.r(rcode))
        rcode = """dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)"""
        print(robjects.r(rcode))
        # add to the path
        rcode = """.libPaths(Sys.getenv("R_LIBS_USER"))"""
        print(robjects.r(rcode))
        utils.install_packages(StrVector(names_to_install))


if __name__ == "__main__":
    test_r()
