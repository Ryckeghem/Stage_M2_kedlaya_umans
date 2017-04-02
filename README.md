# Stage_M2_kedlaya_umans
Implémentation de l'algorithme de composition modulaire de polynôme de 
Kedlaya-Umans, dans le cadre de mon stage de fin de master à l'ANSSI.

L'algorithme de Kedlaya-Umans est quasi-linéaire en le degré des polynômes.
L'implémentation montre que cet algorithme n'est pas viable en pratique, 
principalement car il consomme énormément de mémoire, et car la constante dans 
la complexité en terme d'opération élémentaire est énorme.

L'implémentation est faite en C avec la bibliothèque FLINT.
Les commentaires sont en français.

Mon implémentation dispose d'une FFT dans les corps premiers (Fp) optimisée, et 
d'un type de polynôme multivarié dédié à la composition modulaire. Ce type 
est reconnaissable à l'ajout de "_multi" dans le nom des types standards de 
FLINT.

Mon arborescence pour le code C respecte celle de la bibliothèque FLINT. Mais 
plutôt que de modifier la bibliothèque directement, j'ai préféré créé des 
fichiers qui contiennent les ajouts de ma part. Cela explique que mes fichiers 
pour tracer des courbes ou pour faire des tests de validité contiennent des 
inclusions de fichier ".c".

J'ai laissé le code pour le type fq_nmod_poly qui utilise 
fq_nmod_multi_poly, mais l'implémentation n'est pas valide (pour certains cas) 
uniquement pour le fichier multimodular.c de fq_nmod_multi_poly. Cela est dû au 
fait que rédiger un code spécifique pour ce type de polynôme n'a aucune utilité,
 la chose qui semble la plus intelligente à faire est de créer un fichier 
multimodulaire.c qui caste le type fq_nmod_multi_poly en fq_multi_poly, puis 
d'appeler multimodulaire.c pour le type fq_multi_poly. Il convient ensuite de 
caster le résultat en vecteur de fq_nmod.

