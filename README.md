```
        _           _________ 	
  _____(_)_________/ ____/   | 	
 / ___/ / ___/ ___/ /_  / /| | 	
/ /__/ / /  / /__/ __/ / ___ | 	
\___/_/_/   \___/_/   /_/  |_| 	
                              
```
## Project Information

Team number:xohw-400  	<br />
Project name:circFA		<br />
Date:30/06/2018			<br />
Version of uploaded archive:1	<br />
													<br />
University name: Politecnico di Milano				<br />
Supervisor name: Marco Domenico Santambrogio		<br />
Supervisor e-mail: marco.santambrogio@polimi.it		<br />
Participant: Alberto Zeni							<br />	
Email: alberto.zeni@mail.polimi.it					<br />
Participant: Francesco Peverelli					<br />
Email: francesco1.peverelli@mail.polimi.it			<br />
Participant: Enrico Cabri							<br />
Email: enrico.cabri@mail.polimi.it					<br />
													<br />
Board used: xilinx:aws-vu9p-f1:4ddr-xpr-2pr:4.0		<br />
SDAccel Version: 2017.1.op							<br />

## Project Description	
Brief description of project:						<br />
CircFA is an hardware acceleration of a circular RNA aligner. Circular RNAs are a widespread type of RNA in our genome			
that have been recently discovered to be correlated with many types of carcinogenesis and central neural system pathologies.	
Their identfication on patients can be therefore very useful for specific diagnosis. Our solution implements the most intensive
task of the circular RNA identification on genome, which is the alignment process, on FPGA. The proposed project achieves a speedup of	
1,46x for the alignment process over software.

## Project organization
The binary file of the kernel exceeded the maximum dimension of the zip file (>100MB)	
therefore it has been uploaded on google drive, the link is avilable below.				<br />
The doc directory contains the documentation files.										<br />
The file kernelAdaptive.cpp is the source file for the kernel.							<br />
The file maincl.cpp contains the source file for the host.								<br />
The file kseq.h is a library necessary to compile the host.								<br />
The host_circFA file is the binary of the host file										<br />

Instructions to build and test project				<br />
Step 1:	
Go in the project directory and type				<br />
```
source /xilinx/software/SDx/2017.1.op/settings64.sh
```
to source SDAccel									<br />
Step 2:	<br />	

## Link to the Binary File

https://goo.gl/fk3jYv

## Link Link to YouTube Video
