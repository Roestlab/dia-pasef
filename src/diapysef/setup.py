from setuptools import setup

setup(name='diapysef',
      version='0.1',
      description='Analysis and visualization of DIA PASEF data',
      author='Max Frank, Hannes Roest',
      author_email='max.frank@mail.utoronto.ca',
      license='MIT',
      packages=['diapysef'],
      install_requires=[
          'pandas',
          'numpy',
          'matplotlib',
          'pyopenms'],
      package_data={
          'diapysef': ['data/*']
      },
      scripts=['scripts/get_dia_windows.py',
               'scripts/annotate_mq_ionmobility.py',
               'scripts/plot_dia_windows.py',
               'scripts/convertTDFtoMzML.py'],
      zip_safe=False)
