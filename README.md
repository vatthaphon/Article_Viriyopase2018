# Cooperation and competition of gamma rhythms in delayed pulse-coupled oscillators

<p align="justify">Illustration of interactions between ING and PING using a simple network of delayed pulse-coupled oscillators</p>

<p align="justify">This repository contains the source code accompanying the paper titled "Analyzing the competition of gamma rhythms with delayed pulse-coupled oscillators in phase representation", publised at Physical Review E. 
The goal of this study is to theoretically investigate how ING and PING interact in a simple network of two oscillators using phase representaiton.</p>

## Table of Contents  
[1. Introduction](#Introduction)  
[2. Installation](#Installation)  
[3. Usage](#Usage)  
[4. File Structure](#FileStructure)  
[5. License](#License)  
          
## Introduction<a name="Introduction"/>
<p align="justify">This repository provides the implementation of the methods and algorithms described in the paper "Analyzing the competition of gamma rhythms with delayed pulse-coupled oscillators in phase representation". 
The main goal of this study is to theoretically investigate how ING and PING interact in a simple network of two oscillators. We generally find that we observe the same mechanisms of coperation and competition of gamma rhythms
as we found in larger complex networks in the simpler network.</p>

You can find the study using the larger complex networks [here](https://github.com/vatthaphon/Article_Viriyopase2016).

## Installation<a name="Installation"/>
- Matlab 2018.

## Usage<a name="Usage"/>
- Run the Matlab files in each directory to generate dynamics of the network as varying various parameters.

## File Structure<a name="FileStructure"/>
```plaintext
├── Types11BifDiagFullVaryPhiE/            # Generate dynamics as varying the intrinsic frequency of E oscillator.
├── Types11BifDiagFullVaryPhiEPhiI/        # Generate dynamics as varying both intrinsic frequencies.  
├── Types11BifDiagFullVaryPhiI/            # Generate dynamics as varying the intrinsic frequency of I oscillator.
├── Types11BifDiagFullVarygI2E/            # Generate dynamics as varying the I->E strength.
├── Types11BifDiagFullVarygI2I/            # Generate dynamics as varying the self inhibition strength.
├── Types12BifDiagFullVaryPhiE/            # Generate dynamics as varying the intrinsic frequency of E oscillator.
├── Types12BifDiagFullVaryPhiEPhiI/        # Generate dynamics as varying both intrinsic frequencies.  
├── Types12BifDiagFullVaryPhiI/            # Generate dynamics as varying the intrinsic frequency of I oscillator.
├── Types12BifDiagFullVarygI2E/            # Generate dynamics as varying the I->E strength.
├── Types12BifDiagFullVarygI2I/            # Generate dynamics as varying the self inhibition strength.
├── TypesQIF11BifDiagFullVaryPhiEPhiI/     # Generate dynamics as varying both intrinsic frequencies.  
├── Viriyopase2018.pdf                    # The article.   
├── README.md  
└── LICENSE
```
## License<a name="License"/>
This project is licensed under the MIT License - see the LICENSE file for details.

