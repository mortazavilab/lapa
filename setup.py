from setuptools import setup, find_packages


requirements = [
    'setuptools',
    'tqdm',
    'click',
    'pandas',
    'pybigwig',
    'scipy',
    'bamread>=0.0.9',
    'pyranges>=0.0.71',
    'kipoiseq>=0.3.0',
    'matplotlib'
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest']

setup(
    author="M. Hasan Çelik",
    author_email='muhammedhasancelik@gmail.com',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    description="Alternative polyadenylation detection from diverse data sources",
    install_requires=requirements,
    license="MIT license",
    entry_points='''
        [console_scripts]
        lapa=lapa.main:cli
    ''',
    keywords=['genomics', 'long read RNA-seq', 'APA'],
    name='lapa',
    packages=find_packages(include=['lapa']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mortazavilab/lapa',
    version='0.0.1',
)
