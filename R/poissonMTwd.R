#############################################################
#	poissonMTsetwd, poissonMTgetwd functions
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 18 October 2018
#	Version: 0.1
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################

poissonMTsetwd <- function(path) {
  options("poissonMT:wd" = path)
}

poissonMTgetwd <- function() {
  getOption("poissonMT:wd", NULL)
}
