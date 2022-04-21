""" Primer3 and isPCR installer """

import logging
import os
import shutil
import sys
import subprocess
import tempfile

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
            os.chdir(cwd)
            return output.returncode
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
    # want to do the build in a tmpdir and move executable primer_tk dirXSD
    thisdir = os.path.dirname(os.path.realpath(__file__))
    logging.info(f' checking for presence of isPCR, cloning if not here.')
    homedir = os.path.expanduser('~')
    os.environ['MACHTYPE'] = os.uname().machine
    if not os.path.exists(f"{homedir}/bin/{os.environ['MACHTYPE']}"):
        os.makedirs(f"{homedir}/bin/{os.environ['MACHTYPE']}")
    if not os.path.exists(f'{thisdir}/isPcr33'):
        os.makedirs(f'{thisdir}/isPcr33')
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        logging.info(f" downloading isPCR into {tmpdir}.")
        logging.info(f" cmd: wget --no-check-certificate -O {tmpdir}/isPcr33.zip https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip")
        output = subprocess.run(["wget", "--no-check-certificate", "-O", f'{tmpdir}/isPcr33.zip', "https://hgwdev.gi.ucsc.edu/~kent/src/isPcr33.zip"])
        if output.returncode != 0:
            logging.error(" wget of isPCR failed, exiting.")
            return output.returncode
        logging.info(" wget succeeded, unzipping and removing zip")
        output = subprocess.run(["unzip", f"{tmpdir}/isPcr33.zip", "-d", f'{tmpdir}/isPcrSrc'])
        if output.returncode != 0:
            logging.error(" unzip of isPCR failed. Exiting.")
            return output.returncode
        logging.info(f" building isPCR.")
        try:
            os.chdir(f"{tmpdir}/isPcrSrc")
            os.makedirs(f"{tmpdir}/isPcrSrc/lib/{os.environ['MACHTYPE']}")
            os.chdir(f"{tmpdir}/isPcrSrc/lib")
            logging.info(f" inside {os.getcwd()}")
            output = subprocess.run(["make", 'HG_WARN=""'])
            if output.returncode != 0:
                logging.error(f" make of {tmpdir}/isPcrSrc/lib failed. Exiting")
                return output.returncode
            os.chdir(f"{tmpdir}/isPcrSrc")
            output = subprocess.run(["make", 'HG_WARN=""'])
            if output.returncode != 0:
                logging.error(f" make of {tmpdir}/isPcrSrc failed. Exiting")
                return output.returncode
            shutil.move(f"{homedir}/bin/{os.environ['MACHTYPE']}/isPcr", f"{thisdir}/isPcr33")
            logging.info(' build successful! restoring path.')
            os.chdir(thisdir)
            return 0
        except FileNotFoundError as e:
            logging.error(f" {e}. File was missing. Recovering and exiting")
            os.chdir(cwd)
            return 1
        except PermissionError as e:
            logging.error(f' {e}. you do not have correct permissions to build at {tmpdir}.')
            os.chdir(cwd)
            return 2
        finally:
            logging.info(f' restoring {cwd}')
            os.chdir(cwd)
            return 0

        

def main():
    rc = install_primer3()
    if rc != 0:
        return rc
    rc = install_ispcr()
    return rc

if __name__ == "__main__":
    sys.exit(main())
