#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E.Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

# This function identifies poorly fitted measurements for specific experimental conditions.
# It returns a list of possible indices and mse's pointing to possible connections to be added
# during the feeding process

# Inputs:
# Mandatory:  A cnolist object containing the data (cnolist)
#             A model to optimize (model)
#             A simulation data object for a specific set of dynamic parameters as returned by the getLBodeSimFunction.R function
#
# Optional:   A percentage of measurements identified as misfits based on the MI value (percMI = 0.1 or 10% by default).
#             A number of equally spced points to split the time-course data and observation (nSplines = 100 by default)
#             A binning number (nBins = 10 by default)
#             A spline interpolation method (spliningMethod = "fmm" by default. For more check the spline() function)
#             A thrreshold parameter for minimal misfit to be considered (mseThresh = 0.05 by default)
#             A method about how to identify the misfits: It can be 'mse' or 'mi' (method = 'mse' by default)

identifyMisfitIndices <- function(cnolist = cnolist, model = model, simData = simData, mseThresh = 0.05){
  
  
  indices = computeMSE(cnolist = cnolist, model = model, mseThresh = mseThresh, simData = simData)
  
  return(indices)
  
}