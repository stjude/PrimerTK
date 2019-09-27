## Introduction to PrimerTK

The Primer ToolKit (PrimerTK) seeks to alleviate several of the challenges of primer design such as thermodynamic considerations, avoiding repeats, GC content, and length using a simple, fast, and powerful workflow that will enable users to generate high quality primers for large sets of samples. The workflow for primer generation follows the following design:

```

          ########################                      ####################                     ############################
          # input positions file #                      # reference genome #                     # thermodynamic parameters #
          ########################                      ####################                     ############################
                    | |                                          | |                                          | |
                     |                                            |                                            |
                         |                                        |                                        |
                             |                                    |                                    |
                                 |                                |                                |
                                     |                            |                            |
                                         |                        |                        |
                                         ###################################################
                                         # Primer3 Primer Generation and Scoring + Ranking #
                                         ###################################################
                                                                  |
                                                                  |
                                                                  |
                                                                  |
                                         ###################################################
                                         #        Primer3 output parsing + PCR Setup       #
                                         #     Optional Multiplexing and filtering primers #
                                         ###################################################
                                                                  |
                                                                  |
                                                                  |
                                                                  |
                                         ###################################################
                                         #  In Silico PCR off target amplification checks  #
                                         ###################################################
                                                                  |
                                                                  |
                                                                  |
                                                                  |
                                         ###################################################
                                         #    Filter off target products, return primer    #
                                         #    thermodynamics and pcr product results and   #
                                         #      top ranking final primers by position      #
                                         ###################################################
                                         |                        |                        |
                                     |                            |                            |
                                 |                                |                                |
                        #################                         |                          #####################
                        #  Top Primers  #                         |                          # Plate Order Sheet #
                        #################                         |                          #####################
                                                                  |
                                                                  |
                                                                  |
                                         ###################################################
                                         #   Annotate all primers with SNPs using tabix    #
                                         ###################################################
```

This is a general design schematic of the workflow, more in depth discussion will be covered in each walkthrough:

* [standard and multiplex tutorial](standard_and_multiplex.md)
* [structural variants](structural_variants.md)

Next, proceed to [installation](installation.md).