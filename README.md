# Descriptor_Optimization
This fold contains all necessary code files to optimize the NMBE with the cost functions: Silhoutte, Davies-Bouldin index or Calinski-Harabasz index.

The matching process includes the following four steps:

cost_functions.py.......: calculating the the cost functions, that was clustering validation index in this case
descriptors.py..........: calculating Normalized Multiscale Bending Energy (NMBE)
nmbe_multi.py...........: use descriptors.py script to describe all dataset using multiple process
main.py.................: main script that start the optimization

There are some parameters in the codes, and they should be changed by modifying the codes in specific places.


For more details, please refer to
CARNEIRO, A. C.; LOPES, J. G. F.; ROCHA NETO, J. F. S.; SOUZA, M. M. S.; MEDEIROS,
F. N. S.; BEZERRA, F. N. On the evaluation of cost functions for parameter optimization of a
multiscale shape descriptor. In: IEEE. International Symposium on Signal Processing and
Information Technology (ISSPIT), 2017.


Allan Cordeiro Carneiro
Ms, Instituto Federal de Educação, Ciência e Tecnologia do Ceará
Ceará, Brazil
May 4th, 2018
