from setuptools import setup

setup(
        name='smk_wgbs',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=['smk_wgbs'],
        package_data={
            '': ['*.R', '*.snakefile', '*.yml', '*.yaml', '*.sh', '*.smk', '*.rules'],
        },
        entry_points={
            'console_scripts': [
                'smk_wgbs_snakefile = smk_wgbs.tools:print_snakefile_fp',
                'smk_wgbs_demo_config = smk_wgbs.tools:print_demo_config_fp',
            ]
},
install_requires=[],
        python_requires='>=3.6',
)
