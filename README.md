# indistinguishability-fourier-space
## Setup
* clone this repository to /path/to/repo
* add /path/to/repo  to your $PATH
* add /path/to/repo to your Wolfram init.m file
	This file can be located (normally ~/.Wolfram/Kernel/init.m) using
  ```bash
	wolframscript -code 'Print[FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}]]'
  ```
 	and the path can be added by appending the line
	```
	AppendTo[$Path, "/path/to/repo"];
