import setuptools

setuptools.setup(
    name='primer_tk',
    version='0.2.2',
    scripts=['./scripts/primer_tk'],
    author="Dennis Kennetz",
    author_email="dennis.kennetz@stjude.org",
    description="A toolkit to design primers in multiplex pools and around SVs.",
    packages=['primer_tk', 'primer_tk.tests'],
    package_dir={'primer_tk': 'src/primer_tk', 'primer_tk.tests': 'test/python_tests'},
    package_data={'test': [
        'data/*']},
    install_requires=[
        'setuptools',
        'pandas >= 0.22.0',
        'numpy >= 1.16.0',
        'biopython >= 1.70',
        'pysam==0.15.2'
    ],
    python_requires='>=3.5.*',
    test_suite='test',
    tests_require=['unittest', 'coverage'],
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
