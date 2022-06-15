from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()


requirements = [
    'setuptools',
    'tqdm',
    'click',
    'pandas',
    'pybigwig',
    'scipy',
    'bamread>=0.0.10',
    'pyranges>=0.0.71',
    'sorted-nearest==0.0.33',
    'kipoiseq>=0.3.0',
    'matplotlib',
    'betabinomial',
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest']

setup(
    name='lapa',
    version='0.0.4',

    author="M. Hasan Çelik",
    author_email='muhammedhasancelik@gmail.com',
    url='https://github.com/mortazavilab/lapa',

    keywords=['genomics', 'long read RNA-seq', 'APA'],
    description="Tools for alternative polyadenylation detection and analysis "
    "from diverse data sources (3'seq and long-reads) and "
    "transcript start site detection and analysis from long-read RNA-seq.",

    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    license="MIT license",
    long_description=readme + '\n',
    long_description_content_type='text/markdown',

    install_requires=requirements,
    setup_requires=setup_requirements,

    entry_points='''
        [console_scripts]
        lapa=lapa.main:cli_lapa
        lapa_tss=lapa.main:cli_lapa_tss
        lapa_link_tss_to_tes=lapa.main:cli_lapa_link_tss_to_tes
        lapa_correct_talon_gtf=lapa.main:cli_lapa_correct_talon_gtf
        lapa_correct_talon=lapa.main:cli_lapa_correct_talon
    ''',
    packages=find_packages(include=['lapa*']),
    include_package_data=True,

    test_suite='tests',
    tests_require=test_requirements,
)
