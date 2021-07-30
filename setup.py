from setuptools import setup, find_packages


requirements = [
    'setuptools',
    'tqdm',
    'click',
    'pyranges',
    'pandas',
    'kipoiseq',
    'pybigwig',
    'scipy'
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
    url='https://github.com/MuhammedHasan/longread_postprocessing',
    version='1.0.0',
)
