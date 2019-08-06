import setuptools

setuptools.setup(
    name='primer_tk',
    version='0.2.2',
    scripts=['./scripts/primer_tk'],
    author="Dennis Kennetz",
    author_email="dennis.kennetz@stjude.org",
    description="A toolkit to design primers in multiplex pools and around SVs.",
    packages=['primer_tk', 'genome_standard', 'genome_sv', 'cwl'],
    package_dir={'primer_tk': 'lib/primer_tk',
                 'genome_standard': 'test/genomeStandard',
                 'genome_sv': 'test/genomeSV',
                 'cwl': 'cwl'},
    package_data={'genome_standard': [
        'python_tests/*.py',
        'data/*',
        'bats_tests/commandlinetool/*',
        'bats_tests/data/clt/*',
        'bats_tests/data/wf/*',
        'bats_tests/inputs/clt/*',
        'bats_tests/inputs/wf/*',
        'bats_tests/workflow/*'],
                  'genome_sv': [
                      'python_tests/*.py',
                      'data/*',
                      'bats_tests/commandlinetool/*',
                      'bats_tests/data/sv/*',
                      'bats_tests/data/sv/wf/*',
                      'bats_tests/inputs/sv/*',
                      'bats_tests/inputs/sv/wf/*',
                      'bats_tests/workflow/*'],
                  'cwl': [
                      './*',
                      'tools/*.cwl']},
    include_package_data=True,
    install_requires=[
        'setuptools',
        'pandas >= 0.22.0',
        'numpy >= 1.16.0',
        'biopython >= 1.70',
        'pysam==0.15.2'
    ],
    python_requires='>=3.5.*',
    zip_safe=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        'Intended Audience :: Biology/Genomics',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: OS Independent'
    ],
)
