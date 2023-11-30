#Please contact carolynkoehn@u.boisestate.edu if you have any questions. 

#############################
# Function for running optimiaztions and generating customized output tables
# Author: Carolyn Koehn
# Modified: 2023-02-13
# Function: generate_optim()
# Inputs:
  # inputs: A RasterStack of ecosystem services maps (units in ES per pixel). Optimized solutions will be generated that target a specified amount of these ES.
  # cost: A RasterLayer of cost per pixel. Fed directly into optimization problem with no modification
  # targets: A vector (or single number) specifying protection target(s) to be applied to the ES layers. A vector will cause multiple optimizations to be run. All ES will have the same target applied to them. A matrix may be used where each row is an optimization run and each column is an ES from the inputs RasterStack. This will be used in add_relative_targets() in prioritizr
  # mask: An optional RasterLayer. If table.columns includes "percent_mask_pixels_in_solution", the percent of pixels in the solution that are within the mask layer will be calculated.
  # locked.in and locked.out: see prioritizr help file for accepted inputs. These will be fed directly into the optimization problem with no modifications
  # output.table: If TRUE, a table with selected columns will be returned. If FALSE, no table is generated or returned.
  # output.rasters: IF TRUE, a list with a raster of each optimization solution will be returned. If FALSE, solution rasters will not be saved or returned.
  # table.columns: A vector containing the types of columns that will be included in the output.table.
    # "total.cost": The total cost of the solution in the units of the cost input
    # "percent_ES_protected": The percent of the total ES of each input map that is retained in the optimized solution. Ideally, these would match or exceed the specified target for that run
    # "sum_ES_protected": The total ES protected in the solution. Units are in the units of the input map
    # "total_pixels": Number of pixels included in solution. Can be used to convert to total land area included in solution.
    # "percent_mask_pixels_in_solution": If mask is provided, the percent of solution pixels that are within the mask. Can be used to assess how much land is protected in an area of interest

# Outputs:
  # If output.table=TRUE and output.rasters=TRUE, a list of [1] the output table and [2] a list of solution rasters
  # If one of these options is FALSE, only the option listed as TRUE is provided
  # If both are FALSE, nothing is returned
###################################

generate_optim <- function(inputs, cost, targets, mask=NULL, 
                           locked.in=NULL, locked.out=NULL,
                           output.table=TRUE, output.rasters=FALSE, 
                           table.columns=c("total_cost","percent_ES_protected","sum_ES_protected","total_pixels","percent_mask_pixels_in_solution")) {
  
  # Checks and errors
  if(("percent_mask_pixels_in_solution" %in% table.columns) & is.null(mask) & output.table==TRUE) {
    stop("Must provide mask layer if percent_mask_pixels_in_solution is included in table output")
  }
  
  # Pre-process inputs
  
  #prioritizr cannot handle negative inputs or inputs that are too large. here I modify the raster layers to meet standards and save which transformations I performed
  changes <- data.frame(layer = 1:nlyr(inputs),
                        negative = FALSE,
                        divided_by = 1)
  for(i in 1:nlyr(inputs)) {
    #change negative to 0 and send warning message
    if(min(inputs[[i]][],na.rm=TRUE) < 0) {
      inputs[[i]][inputs[[i]] < 0] <- 0
      changes$negative[i] <- TRUE
      warning(paste("Negative values in ES layer",i,"converted to 0"))
    }
    #divide layers so maximum is <1000 and keep track of which ones
    if(max(inputs[[i]][],na.rm=TRUE) > 1000) {
      digits_before_decimal <- floor(log10(max(inputs[[i]][],na.rm=TRUE))) + 1
      divide_by <- 1 * 10^(digits_before_decimal - 2)
      inputs[[i]] <- inputs[[i]]/divide_by
      changes$divided_by[i] <- divide_by
    }
  }
  
  print("Pre-processing complete")
  
  # Create containers for optimization outputs
  
  if(output.table == TRUE) {
    if(is.matrix(targets)) {
      if(ncol(targets) != nlyr(inputs)) {
        stop("Number of columns in targets matrix must equal number of input layers")
      }
      target_columns <- paste("target_ES",1:nlyr(inputs),sep="")
      out_table <- data.frame(targets)
      colnames(out_table) <- target_columns
    } else {
      out_table <- data.frame(target=targets)
    }
    
    if("total_cost" %in% table.columns) {
      out_table <- cbind(out_table, total_cost=NA)
    }
    if("percent_ES_protected" %in% table.columns) {
      column_names1 <- paste("percent_ES",1:nlyr(inputs),"_protected",sep="")
      add <- matrix(data=NA, nrow=nrow(out_table), ncol=length(column_names1))
      colnames(add) <- column_names1
      out_table <- cbind(out_table, add)
    }
    if("sum_ES_protected" %in% table.columns) {
      column_names2 <- paste("sum_ES",1:nlyr(inputs),"_protected",sep="")
      add <- matrix(data=NA, nrow=nrow(out_table), ncol=length(column_names2))
      colnames(add) <- column_names2
      out_table <- cbind(out_table, add)
    }
    if("total_pixels" %in% table.columns) {
      out_table <- cbind(out_table, total_pixels=NA)
    }
    if("percent_mask_pixels_in_solution" %in% table.columns) {
      out_table <- cbind(out_table, percent_mask_pixels_in_solution=NA)
    }
  }
  
  if(output.rasters == TRUE) {
    out_rasts <- vector("list", NROW(targets))
  }
  
  # Set up optimization problems loop
  
  for(i in 1:NROW(targets)) { #generate solutions for the range of protection targets
    problem <- problem(cost,inputs) %>% #create problem with cost and input layers
      add_min_set_objective() %>% #find solution with lowest cost
      add_relative_targets(targets[i]) %>% #set yield protection targets
      add_gurobi_solver(verbose=FALSE)
    
    if(!is.null(locked.in)) {
      problem <- problem %>%
        add_locked_in_constraints(locked.in)
    }
    
    if(!is.null(locked.out)) {
      problem <- problem %>%
        add_locked_out_constraints(locked.out)
    }
    
    solution <- solve(problem) #use Gurobi to find optimal solution
    
    # fill output objects
    
    if(output.table == TRUE) {
      if("total_cost" %in% table.columns) {
        out_table$total_cost[i] <- as.numeric(eval_cost_summary(problem,solution)["cost"])
      }
      if("percent_ES_protected" %in% table.columns) {
        out_table[i,column_names1] <- t(eval_feature_representation_summary(problem,solution)[,'relative_held'])*100
      }
      if("sum_ES_protected" %in% table.columns) {
        out_table[i,column_names2] <- t(eval_feature_representation_summary(problem,solution)[,'absolute_held'])
        # convert back to original scale
        out_table[i,column_names2] <- out_table[i,column_names2] * changes$divided_by
      }
      if("total_pixels" %in% table.columns) {
        out_table$total_pixels[i] <- as.numeric(eval_n_summary(problem,solution)["n"])
      }
      if("percent_mask_pixels_in_solution" %in% table.columns) {
        out_table$percent_mask_pixels_in_solution[i] <- as.numeric(table(solution[!is.na(mask)])[2]/subset(freq(mask), freq(mask)[,"value"]==1)[,"count"])*100
      }
    }
    
    if(output.rasters==TRUE) {
      out_rasts[[i]] <- solution
    }
    print(paste("Optimization ",i,"/",NROW(targets)," complete",sep=""))
  }
  
  # return requested output
  if(output.table==TRUE & output.rasters==TRUE) {
    return(list(out_table, out_rasts))
  } else if(output.table==TRUE & output.rasters==FALSE) {
    return(out_table)
  } else if(output.rasters==TRUE & output.table==FALSE) {
    return(out_rasts)
  }
}

