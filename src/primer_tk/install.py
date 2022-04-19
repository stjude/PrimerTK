""" Primer3 and isPCR installer """

import logging
import os
import sys
import subprocess

def install_primer3():
    """ Installs primer3 in the program, so it can be used directly by the program """
    thisdir = os.path.dirname(os.path.realpath(__file__))
    logging.info(f' checking for presence of primer3, cloning if not here.')
    if not os.path.exists(f'{thisdir}/primer3_install'):
        output = subprocess.run(["git", "clone", "https://github.com/primer3-org/primer3.git", f'{thisdir}/primer3_install'])
        if output.returncode != 0:
            logging.error(' clone of primer3 failed. check git is in path and you can clone')
            return output.returncode
    elif os.path.exists(f'{thisdir}/primer3_install/src/primer3_core'):
        logging.warning(f' primer3 appears to already be built at {thisdir}/primer3_install/src/primer3_core. Nothing to do!')
        return 0
    elif os.path.exists(f'{thisdir}/primer3_install') and not os.path.exists(f'{thisdir}/primer3_install/src/primer3_core'):
        logging.warning(f' primer3 appears to have been cloned but not built. Building now!')
    else:
        logging.error(f' something went wrong!')
        return 1
    cwd = os.getcwd()
    try:
        logging.info(' building primer3')
        os.chdir(f'{thisdir}/primer3_install/src')
        output = subprocess.run(["make"])
        if output.returncode != 0:
            logging.error(' build of primer3 failed. Exiting')
            sys.exit(output)
        logging.info(' build successful! restoring path.')
        os.chdir(thisdir)
    except FileNotFoundError:
        logging.error(f' {thisdir}/primer3_install/src did not exist as a path.')
        os.chdir(cwd)
        return 1
    except PermissionError:
        logging.error(f' you do not have correct permissions to build at {thisdir}/primer3_install/src.')
        os.chdir(cwd)
        return 2
    finally:
        logging.info(f' restoring {cwd}')
        os.chdir(cwd)
        return 0

def install_ispcr():
    """ Installs ispcr directly into the program, so it can be used by the program """
    # want to do the build in a tmpdir and move executable primer_tk dir
    thisdir = os.path.dirname(os.path.realpath(__file__))
    logging.info(f' checking for presence of isPCR, cloning if not here.')
    if not os.path.exists(f'{thisdir}/isPcrSrc'):
        logging.info(" getting isPCR.")
        logging.info(" cmd: wget --no-check-certificate https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip")
        output = subprocess.run(["wget", "--no-check-certificate", "-O", f'{thisdir}/isPcr33.zip', "https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip"])
        if output.returncode != 0:
            logging.error(" wget of isPCR failed, exiting.")
            return output.returncode
        logging.info(" wget succeeded, unzipping and removing zip")
        output = subprocess.run(["unzip", f"{thisdir}/isPcr33.zip", "-d", f'{thisdir}/isPcrSrc'])
        if output.returncode != 0:
            logging.error(" unzip of isPCR failed. Exiting.")
            return output.returncode
        

def main():
    rc = install_primer3()
    if rc != 0:
        return rc
    rc = install_ispcr()
    return rc

if __name__ == "__main__":
    sys.exit(main())
