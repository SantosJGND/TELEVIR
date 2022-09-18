<p align="center"><img src="logo/logo_insaflu.png" alt="INSaFLU" width="300"></p>


[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)


# Televir - docker installation
## Hardware Requirements

* Processor: 8 cores (4 minimal);
* RAM: 32GB of memory (16GB minimal);
* Disk Space: 512GB (suggestion; depends on the volume of data to process);


## Installation

Docker:

* Install [docker](https://docs.docker.com/engine/install/) in your linux server;
* Install [docker-compose](https://docs.docker.com/compose/install/) in your linux server;

	$ sudo curl -L "https://github.com/docker/compose/releases/download/1.27.4/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
	
	$ sudo chmod +x /usr/local/bin/docker-compose

* Install docker extensions [local-persist](https://github.com/MatchbookLab/local-persist);

	$ curl -fsSL https://raw.githubusercontent.com/MatchbookLab/local-persist/master/scripts/install.sh > install.sh
	
	$ chmod a+x install.sh

	$ sudo ./install.sh

:warning: If you're uncomfortable running a script you downloaded off the internet with sudo, you can extract any of the steps out of the install.sh script and run them manually.

INSaFLU:

	$ git clone https://github.com/INSaFLU/docker.git
	$ cd docker
	
	## to define the directory where the data will be saved and the web port exposed, edit the .env file: 
	$ cp .env_temp .env
	$ vi .env
	OR
	$ nano .env
	
	## add your user account to docker group to use docker without sudo
	$ sudo usermod -aG docker $USER
	$ sudo chmod 666 /var/run/docker.sock
	
	## test if everything is OK
	$ docker ps
	$ docker run hello-world 
	
	## build INSaFLU
	$ ./build.sh
	$ ./up.sh
	
	## create an user, in other terminal or you can use 'screen' in previous steps
	$ docker exec -it insaflu-server create-user
	
Now, you can go to a web explorer and link to the address "127.0.0.1:<port defined in .env>". Default port is 8080

To stop:

	$ ./stop.sh

To start again:

	$ ./up.sh
	
