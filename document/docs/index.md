# Welcome to FederEI

FederEI is a federated library matching solution for EI-MS-based compound identification. 
It establishes a server-to-server connection framework seamlessly integrated into a user-friendly front-end software. 
By keeping data localized within each laboratory's server, FederEI minimizes the need for sharing 
sensitive spectral information across multiple entities, thus mitigating privacy concerns.    

Here is the workflow: The user submits mass spectrometry data through the user interface (UI), 
and FederEI dispatches it to the central-server. Upon reception, the central-server distributes 
the file to all client-server(s) based on their respective IP addresses stored in the client-server registry. 
Subsequently, each client-server searches for results in its local database using library matching algorithm 
and transmits them back to the central-server. The central-server tallies the received search results for the 
file until all client-servers have responded. Finally, central-server organizes and summarizes the results, 
and return the final results to the corresponding user.

[![fig1.jpg](https://i.postimg.cc/ZYs9Pc2C/fig1.jpg)](https://postimg.cc/WtrbTM2v)

## Documentation

- [Installation and Deployment](https://hcji.github.io/FederEI/Installation_and_Deployment/)
- [User Guide](https://hcji.github.io/FederEI/User_Guide/)
- [Advanced Usage](https://hcji.github.io/FederEI/Advanced_Usage/)

## Sources on GitHub
[FederEI GitHub](https://github.com/hcji/FederEI)

## Changelog

- verision 1.0: *first release* at 2024.04.30

## Citations
- [Yang, Q., Ji, H., Xu, Z. et al. Ultra-fast and accurate electron ionization mass spectrum matching for compound identification with million-scale in-silico library. Nat Commun 14, 3722 (2023)](https://doi.org/10.1038/s41467-023-39279-7)


## Contact
- jihongchao@caas.cn
