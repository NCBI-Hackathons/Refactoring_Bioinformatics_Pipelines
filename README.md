Cherry-picked, refactored Nextflow pipelines for speed, picked from
https://github.com/nextflow-io/awesome-nextflow for modification to
scale up.

For information on the domain-specific language used see https://github.com/nextflow-io/nextflow


Analysis of NGS data involves multiple technical steps such as installation of the software components of bioinformatics pipelines; coordinating format conversions and data flow between pipeline components; managing software versions and updates; automating execution for multiple runs; supplying the required computational and data storage infrastructure; and last but not least, providing an intuitive user interface for non-bioinformatics experts. To overcome these challenges, bioinformatics software developers have leveraged technologies such as virtual machines and Docker containers ([1], [2]) for distributing preconfigured bioinformatics software that can run on any computational platform. The use of virtualization technology and cloud computing saves significant development time and cost, as the software does not need to be set up from scratch at each laboratory. The increased interest for applications of virtualization for NGS data analysis is evident through many recent studies, ranging from comparing performance of virtual machines to conventional computing [3], and bioinformatics-specific Docker container repositories [4]. 
 
We propose that by leveraging the high-level, domain-specific language Nextflow [ref] that is designed for bioinformatics, we can develop a novel model called Science as a Service (SciaaS). Following this model, bioinformatics software developers can easily deploy software and data analysis on a variety of computing environments ranging from lab computers, data center clusters setups and cloud computers. With a SciaaS model for bioinformatics, researchers can easily analyze sequencing data, and remove the bottleneck in smaller laboratories, as well as provide easy access for large-scale computational and storage capabilities with the use of cloud computing.


References
1. Krampis K, Booth T, Chapman B, Tiwari B, Bicak M, Field D, et al. Cloud BioLinux: pre-configured and on-demand bioinformatics computing for the genomics community. BMC Bioinformatics 2012;13:1â 8. 
 
2. Hosny A, Vera-Licona P, Laubenbacher R, Favre T. AlgoRun: A Docker-based packaging system for platform-agnostic implemented algorithms. Bioinformatics. 2016;32:2396â8. 
 
3. Di Tommaso P, Palumbo E, Chatzou M, Prieto P, Heuer ML, Notredame C. The impact of Docker containers on the performance of genomic pipelines. PeerJ [Internet]. 2015;3:e1273. 
 
4. Moreews F, Sallou O, MÃ©nager H, Le bras Y, Monjeaud C, Blanchet C, et al. BioShaDock: a community driven bioinformatics shared Docker-based tools registry. F1000Research [Internet]. 2015;1â 9. 

