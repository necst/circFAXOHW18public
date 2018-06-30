```
        _           _________ 	
  _____(_)_________/ ____/   | 	
 / ___/ / ___/ ___/ /_  / /| | 	
/ /__/ / /  / /__/ __/ / ___ | 	
\___/_/_/   \___/_/   /_/  |_| 	
                              
```
## Project Information

Team number: xohw-400  	<br />
Project name: circFA		<br />
Date: 30/06/2018			<br />
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
CircFA is an hardware acceleration of a circular RNA aligner. Circular RNAs are a widespread type of RNA in our genome			
that have been recently discovered to be correlated with many types of carcinogenesis and central neural system pathologies.	
Their identfication on patients can be therefore very useful for specific diagnosis. Our solution implements the most intensive
task of the circular RNA identification on genome, which is the alignment process, on FPGA. The proposed project achieves a speedup of 1,46x for the alignment process over software.

## Project organization
The binary file of the kernel exceeded the maximum dimension of the zip file (>100MB)	
therefore it has been uploaded on google drive, the link is avilable below.				<br />
The doc directory contains the documentation files.										<br />
The file alignment_sw_reference.cpp contains the reference version of the alignment algorithm in software <br />
The file kernelAdaptive.cpp is the source file for the kernel.							<br />
The file maincl.cpp contains the source file for the host.								<br />
The file kseq.h is a library necessary to compile the host.								<br />
The host_circFA file is the binary of the host file.										<br />

Instructions to build and test project				<br />	
Go in the project directory and type in the terminal			<br />
```
source /xilinx/software/SDx/2017.1.op/settings64.sh
```
to source SDAccel	.								<br />
To run  a software emulation type in the terminal:                   
```
make emulation TARGET=sw_emu
```                                                                                     
To run  an hardware  emulation type in the terminal:                   
```
make emulation TARGET=hw_emu
```
To compile the kernel and generate the bitstream type in the terminal:
```
make xclbin TARGET=hw
```
The bitstream will be inside the hw folder.                    <br />
To compile the host type in the terminal:
```
make host TARGET=hw
```
The host willl be inside the hw folder.                             <br />
In order to run the kernel the creation of an F1 instance is needed, 
a complete guide for that can be found here:                  <br />
https://github.com/Xilinx/SDAccel_Examples/wiki/Create,-configure-and-test-an-AWS-F1-instance   <br />
For this project we used an f1.2xlarge instance.              <br />
When both host and the kernel bitstream for AWS have been generated, in order to run the kernel on FPGA type on terminal
```
./host_circFA kernel.awsxclbin ref.fasta query.fasta
```
Where ref.fasta and query.fasta are the two files that will be generated by the host and that will be alligned.
Each files will contain a sequence of 3024 characters that will be alligned in order to test the maximum limit of our application.

## Link to the Binary File
https://goo.gl/fk3jYv

## Link to YouTube Video
https://youtu.be/nHDAhXOTM4A

## Link to Github repository
https://github.com/albertozeni/circFAXOHW18public
