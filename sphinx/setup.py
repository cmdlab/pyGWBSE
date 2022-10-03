from setuptools import setup


setup(
   name='pyGWBSE',
   version='1.0',
   description='python workflow for GW-BSE calculation',
   author='Tathagata Biswas',
   author_email='tbiswas3@asu.edu',
   packages=['pyGWBSE'],
   install_requires=['FireWorks>=1.4.0', 'pymatgen>=2018.6.11',            
                          'monty>=1.0.2', 'atomate>=0.8.4'],   
   package_data={"pyGWBSE": ["inputset.yaml"]}

)
