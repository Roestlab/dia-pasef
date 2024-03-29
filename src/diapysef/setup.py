from setuptools import setup, find_packages

def readme():
    with open('PYPI_README.rst') as f:
        return f.read()

desc = """\
Analysis, conversion and visualization of diaPASEF data."""

setup(name='diapysef',
      version='1.0.10',
      description=desc,
      long_description=readme(),
      long_description_content_type='text/x-rst',
      url="https://github.com/Roestlab/dia-pasef",
      author='Max Frank, Annie Ha, Justin Sing, Josh Charkow, Hannes Roest',
      author_email='hannes.rost@utoronto.ca',
      license='MIT',
      packages=['diapysef'],
      install_requires=[
          'pandas',
          'numpy',
          'matplotlib',
          'seaborn',
          'statsmodels',
          'scipy',
          'scikit-image',
          'pyopenms',
          'patsy',
          'tqdm',
          'joblib',
          'click'],
      package_data={
          'diapysef': ['data/*'],
          '': ['PYPI_README.rst']
      },
      include_package_data=True,
      entry_points={
          'console_scripts': [
              "diapysef=diapysef.main:cli",
              ]
      },
      scripts=['scripts/get_dia_windows.py',
               'scripts/annotate_mq_ionmobility.py',
               'scripts/plot_dia_windows.py',
               'scripts/create_library.py',
               'scripts/high_precision_irt.py'],
      zip_safe=False,
    )
