#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "util.h"
#include "vpr_types.h"
#include "globals.h"
#include "route_export.h"
#include "route_common.h"
#include "route_breadth_first.h"

/********************* Subroutines local to this module *********************/

static boolean breadth_first_route_net(int inet, float bend_cost);

static void breadth_first_expand_trace_segment(struct s_trace *start_ptr, int remaining_connections_to_sink);

static void breadth_first_expand_neighbours(int inode, float pcost, int inet, float bend_cost);

static void breadth_first_add_source_to_heap(int inet);
			
static void multilevel_paritioning_nets(int *, int , int *, int *, int *, int *, int *, int *, int, int); 

static void record_routing_inet_inode_item(struct s_trace *routing_net_start, int *routing_inet_inode);

static int predict_used_inode(struct s_trace *routing_segment_start);

int comparefunc (const void * a, const void * b);

/************************ Subroutine definitions ****************************/

boolean try_breadth_first_route(struct s_router_opts router_opts, t_ivec ** clb_opins_used_locally, int width_fac) {

	/* Iterated maze router ala Pathfinder Negotiated Congestion algorithm,  *
	 * (FPGA 95 p. 111).  Returns TRUE if it can route this FPGA, FALSE if   *
	 * it can't.                                                             */
	float pres_fac;
	boolean success, is_routable, rip_up_local_opins;
	int itry, inet, inode;
	/* Information about routing resource graph */
//	printf("num_nets = %d \n", num_nets);
//  printf("nx = %d ny = %d \n", nx, ny); 
	/* Parallel variable setting in FOUR instances */
	int i = 0, j = 0, k = 0;

	int cross_pid = 0, right_pid = 0;
	int pid, num_pids; //pid is NO. instance, num_pids is the number of instances.
	/* Information about instances required */
   	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pids);		    
//	printf("[Debug %d]------------------num_rr_nodes = %d ---------------- \n", pid, num_rr_nodes);

	MPI_Status status_100, status_101;
	MPI_Status status_50, status_51, status_52;
    MPI_Status status_150, status_151, status_152, status_153;

	int x_median         = (nx+2);
	int y_median         = (ny+2);

	int work_id          = 0;

	/* the two instances to run... 
	int num_work_ids     = 3;
	int owner_pid[]      = {0, 0, 1};
	int first_work_ids[] = {0, 2};
	int dimension[]      = {0};
	int subdimension[]   = {1};
	int median_cutline[]    = {x_median/2};
	int submedian_cutline[] = {y_median/2};
    */	

	/* the four instances to run...  
	int num_work_ids     = 7;
	int owner_pid[]      = {0, 0, 1, 0, 2, 1, 3};
	int first_work_ids[] = {0, 2, 4, 6};
	int dimension[]      = {0, 1, 1};
	int subdimension[]   = {1, 0, 0};
	int median_cutline[]    = {x_median/2, y_median/2, y_median/2};
    int submedian_cutline[] = {y_median/2, x_median/2, x_median/2};    
    */ 

    /* the eight instances to run... 
    int num_work_ids     = 15;
    int owner_pid[]      = {0, 0, 1, 0, 2, 1, 3, 0, 4, 2, 5, 1, 6, 3, 7};
    int first_work_ids[] = {0, 2, 4, 6, 8, 10, 12, 14};
   
   	int dimension[]         = {0, 1, 1, 0, 0, 0, 0};
   	int subdimension[]      = {1, 0, 0, 1, 1, 1, 1};
   	
	int median_cutline[]    = {x_median/2, y_median/2, y_median/2, x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4};
	int submedian_cutline[] = {y_median/2, x_median/2, x_median/2, y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4};
    */

	/* the sixteen instances to run...  
	int num_work_ids     = 31;
	int owner_pid[]      = {0, 0, 1, 0, 2, 1, 3, 0, 4, 2, 5, 1, 6, 3, 7, 0, 8, 4, 9, 2, 10, 5, 11, 1, 12, 6, 13, 3, 14, 7, 15};
	int first_work_ids[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};     // NO. of first work node.

	int dimension[]         = {0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
	int subdimension[]      = {1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
	
	int median_cutline[]    = {x_median/2, y_median/2, y_median/2, x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4, 
	                           y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4, y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4};
   
	int submedian_cutline[]  = {y_median/2, x_median/2, x_median/2, y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4, 
	                           x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4, x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4};
    */ 

    /* the thirty-two instances to run... */	
	int num_work_ids     = 63;
	int owner_pid[]      = {0, 0, 1, 0, 2, 1, 3, 0, 4, 2, 5, 1, 6, 3, 7, 0, 8, 4, 9, 2, 10, 5, 11, 1, 12, 6, 13, 3, 14, 7, 15, 
		                    0, 16, 8, 17, 4, 18, 9, 19, 2, 20, 10, 21, 5, 22, 11, 23, 1, 24, 12, 25, 6, 26, 13, 27, 3, 28, 14, 29, 7, 30, 15, 31}; // store the NO. of instances to num_work_ids.
	int first_work_ids[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62};     // NO. of first work node.
	
	int dimension[]         = {0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // 0:x-coordination 1:y-coordination
	int subdimension[]      = {1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // 0:x-coordination 1:y-coordination
	
	int median_cutline[]    = {x_median/2, 
		                       y_median/2, y_median/2, 
							   x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4, 
		                       y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4, y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4,
		                       x_median/8, x_median/8, (x_median*3)/8, (x_median*3)/8, x_median/8, x_median/8, (x_median*3)/8, (x_median*3)/8, 
							   (x_median*5)/8, (x_median*5)/8, (x_median*7)/8, (x_median*7)/8, (x_median*5)/8, (x_median*5)/8, (x_median*7)/8, (x_median*7)/8}; 
    
	int submedian_cutline[] = {y_median/2, 
		                       x_median/2, x_median/2, 
							   y_median/4, y_median/4, (y_median*3)/4, (y_median*3)/4, 
		                       x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4, x_median/4, x_median/4, (x_median*3)/4, (x_median*3)/4,
		                       y_median/8, y_median/8, (y_median*3)/8, (y_median*3)/8, y_median/8, y_median/8, (y_median*3)/8, (y_median*3)/8, 
							   (y_median*5)/8, (y_median*5)/8, (y_median*7)/8, (y_median*7)/8, (y_median*5)/8, (y_median*5)/8, (y_median*7)/8, (y_median*7)/8};
    
	/****************************** Construct sendrecv variants when parallel router completed overall signals *****************************************************/
	int current_instance_sum_nodes              = 0;
	int current_instance_net                    = 0;

	int current_instance_total_count_nodes      = 0;
	int current_instance_total_count_nets       = 0;
	
	int current_instance_count_net_nodes        = 0;

	int recv_current_instance_sum_nodes         = 0;
	int recv_current_instance_net               = 0;

	int recv_current_instance_total_count_nodes = 0;
	int recv_current_instance_total_count_nets  = 0;
	/***Allocate a new array record pid item to each signal such that information collection when route end***/	
	struct s_trace *new_tptr = 0;
	struct s_trace *recv_new_tptr = 0;
	struct s_route *new_transfer = (t_route *)malloc(sizeof(t_route)*num_rr_nodes);
	assert(new_transfer);
	memset(new_transfer, 0, sizeof(t_route)*num_rr_nodes);

    int isignal = 0;
	int *count_isignal = (int *)malloc(sizeof(int)*num_pids);
    assert(count_isignal);
    memset(count_isignal, 0, sizeof(int)*num_pids);

    int current_instance_total_num_signals = 0;
	int *record_current_instance_each_net_tag = (int *)malloc(sizeof(int)*num_nets); /* updated */
	assert(record_current_instance_each_net_tag);
	memset(record_current_instance_each_net_tag, 0, sizeof(int)*num_nets);

	int *record_current_instance_count_net_nodes = (int *)malloc(sizeof(int)*num_nets);
	assert(record_current_instance_count_net_nodes);
	memset(record_current_instance_count_net_nodes, 0, sizeof(int)*num_nets);
	/** Updated overall rrgraph nodes cost to ensure each instance's cost is same in next iteration if routing false **/
	int *updated_overall_routing_resource_graph_cost = (int *)malloc(sizeof(int)*num_rr_nodes);
	assert(updated_overall_routing_resource_graph_cost);
	memset(updated_overall_routing_resource_graph_cost, 0, sizeof(int)*num_rr_nodes);
    
	/* Allocate the dynamic array in size of num_nets */
	int *cur_nets = 0;
	int *overall_signals = (int *)malloc(sizeof(int)*num_nets);
	assert(overall_signals);
	for (i = 0; i < num_nets; i++) {
		overall_signals[i] = i;
	}

	int cur_count_nets = 0;
	int count_left     = 0;
	int count_right    = 0;
	int count_cross    = 0;
	
	int *left  = (int *)malloc(sizeof(int)*num_nets);
	int *right = (int *)malloc(sizeof(int)*num_nets);
	int *cross = (int *)malloc(sizeof(int)*num_nets);
	
	assert(left);
	memset(left,  0, sizeof(int)*num_nets);
	assert(right);
	memset(right, 0, sizeof(int)*num_nets);
	assert(cross);
	memset(cross, 0, sizeof(int)*num_nets);
   
	/********************************* Here is varieties about the cross nets resegmentation ************************************/
	int count_subleft  = 0;
	int count_subright = 0;
	int count_subcross = 0;

	int *subleft  = (int *)malloc(sizeof(int)*num_nets); 
	int *subright = (int *)malloc(sizeof(int)*num_nets); 
	int *subcross = (int *)malloc(sizeof(int)*num_nets); 
	
    assert(subleft);
	memset(subleft, 0, sizeof(int)*num_nets);
    assert(subright);
	memset(subright, 0, sizeof(int)*num_nets);
    assert(subcross);
	memset(subcross, 0, sizeof(int)*num_nets);

	int count_subcross_subleft_array = 0;
	int count_subright_right_array   = 0;
	int *subcross_subleft_array = (int *)malloc(sizeof(int)*num_nets);
	int *subright_right_array   = (int *)malloc(sizeof(int)*num_nets);

	assert(subcross_subleft_array);
	memset(subcross_subleft_array, 0, sizeof(int)*num_nets);	
	assert(subright_right_array);
	memset(subright_right_array, 0, sizeof(int)*num_nets);

	/* Variable to template paral result. */
	int recv_count_routing_inet_inode = 0, recv_count_right = 0;
	
	/**************Focus inside information send and receive **************************/
    int before_node = 0, after_node = 0, interleaved_node = 0;

	int count_before_routing_inet_inode        = 0;
	int count_after_routing_inet_inode         = 0;
	int count_interleaved_routing_inet_inode   = 0;
	int *record_before_routing_inet_inode      = (int *)malloc(sizeof(int)*num_rr_nodes);
	int *record_after_routing_inet_inode       = (int *)malloc(sizeof(int)*num_rr_nodes);

	int *record_interleaved_routing_inet_inode = (int *)malloc(sizeof(int)*num_rr_nodes);
	int *record_interleaved_delta_inode_occ    = (int *)malloc(sizeof(int)*num_rr_nodes);
    
	int count_routing_inet_inode  = 0;
	int *store_routing_inet_inode = (int *)malloc(sizeof(int)*num_rr_nodes);
	int *store_delta_inet_inode   = (int *)malloc(sizeof(int)*num_rr_nodes);

	assert(record_before_routing_inet_inode); 
	memset(record_before_routing_inet_inode, 0, sizeof(int)*num_rr_nodes);	
	assert(record_after_routing_inet_inode);
	memset(record_after_routing_inet_inode, 0, sizeof(int)*num_rr_nodes);
	assert(record_interleaved_routing_inet_inode);
	memset(record_interleaved_routing_inet_inode, 0, sizeof(int)*num_rr_nodes);
	assert(record_interleaved_delta_inode_occ);
	memset(record_interleaved_delta_inode_occ, 0, sizeof(int)*num_rr_nodes);
	assert(store_routing_inet_inode);
	memset(store_routing_inet_inode, 0, sizeof(int)*num_rr_nodes);
	assert(store_delta_inet_inode);
	memset(store_delta_inet_inode, 0, sizeof(int)*num_rr_nodes);
	/*************** Judge whether instance is quit ***********************************/
	int is_success = 0;
	pres_fac = router_opts.first_iter_pres_fac;
	for (itry = 1; itry <= router_opts.max_router_iterations; itry++) {

		printf("[itry %d][debug %d]-----------------------start-------------------------------\n", itry, pid);
        /**********************************************Receive delta occ and right signals set ******************************************************/
		recv_count_routing_inet_inode = 0;
		recv_count_right = 0;
		if (pid != 0) {
			/****************************Receive (itry-1)th iteration information to routing resource graph ******************************************/
			if (itry > 1) {
				/*********** Judge whether other instance is quit *********************/
				MPI_Bcast(&is_success, 1, MPI_INT, 0, MPI_COMM_WORLD); /* updated */
				if (is_success == 1) {
					printf("[itry %d][Debug %d]NO.%dth instance is quiting...\n", itry, pid, pid);
					MPI_Finalize();
					exit(0);
				}
				/************* Recv routing information from 0th instance *********************************************************/
				
				MPI_Bcast(updated_overall_routing_resource_graph_cost, num_rr_nodes, MPI_INT, 0, MPI_COMM_WORLD); /* updated */
			    for (inode = 0; inode < num_rr_nodes; inode++) {
			    	rr_node[inode].occ = updated_overall_routing_resource_graph_cost[inode]; 
			    }
		    }
			/*****************************************************************************************************************************************/
	//		printf("[itry %d][Debug %d]", itry, pid);
			cross_pid = owner_pid[(first_work_ids[pid]-1)/2];
	//      printf("wait for %d instance\n", cross_pid);
	        /* Receive the occupied node identifier, number and cost. */
			MPI_Recv(&count_routing_inet_inode, 1, MPI_INT, cross_pid, 50, MPI_COMM_WORLD, &status_50);	
			recv_count_routing_inet_inode = count_routing_inet_inode;
			MPI_Recv(store_routing_inet_inode, recv_count_routing_inet_inode, MPI_INT, cross_pid, 51, MPI_COMM_WORLD, &status_51);
            MPI_Recv(store_delta_inet_inode, recv_count_routing_inet_inode, MPI_INT, cross_pid, 52, MPI_COMM_WORLD, &status_52);
			/* Start to update the rr_graph[inode].occ for routing resource graph. */
			for (i = 0; i < recv_count_routing_inet_inode; i++) {
				rr_node[store_routing_inet_inode[i]].occ += store_delta_inet_inode[i]; 
			}

			/* Receive right net and number to assignment scheduling */
			MPI_Recv(&count_subright_right_array, 1, MPI_INT, cross_pid, 100, MPI_COMM_WORLD, &status_100);
			recv_count_right = count_subright_right_array;
			MPI_Recv(subright_right_array, recv_count_right, MPI_INT, cross_pid, 101, MPI_COMM_WORLD, &status_101);
            printf("[itry %d][Debug %d]success receive from %d.\n", itry, pid, cross_pid);
			
			cur_nets = subright_right_array;
			cur_count_nets = recv_count_right;
		    /*
			printf("[itry %d][Debug %d]route right %d nets.\n", itry, pid, cur_count_nets);
		    printf("[itry %d][Debug %d]right net No. is: ", itry, pid);
			for (i = 0; i < cur_count_nets; i++) {
				printf("%d ", cur_nets[i]);
			}
			printf("\n");
			*/
		}
		else {
            /******************************************************************************************************************/
			if (itry > 1) {
				MPI_Bcast(&is_success, 1, MPI_INT, 0, MPI_COMM_WORLD); /* updated */
				if (is_success == 1) {
		            printf("[itry %d][Debug %d]NO.%dth instance is quiting...\n", itry, pid, pid);	
			        /* updated here */		
					free(new_transfer);

					free(overall_signals);
					free(left);
					free(right);
					free(cross);

					free(subleft);
					free(subright);
					free(subcross);
					free(subcross_subleft_array);
					free(subright_right_array);

					free(record_before_routing_inet_inode);      
					free(record_after_routing_inet_inode);       
					free(record_interleaved_routing_inet_inode); 
					free(record_interleaved_delta_inode_occ); 

					free(store_routing_inet_inode);
					free(store_delta_inet_inode);

					free(count_isignal);
					free(record_current_instance_each_net_tag);

					free(updated_overall_routing_resource_graph_cost);
                   
				    return (TRUE);
				    //MPI_Finalize();	
				    //break;
				}
		        /***************** Send rrgraph[inode].occ to next iteration to show same information for all instances ************/
				for (inode = 0; inode < num_rr_nodes; inode++) {
					updated_overall_routing_resource_graph_cost[inode] = rr_node[inode].occ;
				}
				MPI_Bcast(updated_overall_routing_resource_graph_cost, num_rr_nodes, MPI_INT, 0, MPI_COMM_WORLD); /* updated */

            }
            /******************************************************************************************************************/
			cur_nets = overall_signals; /* updated */
			cur_count_nets = num_nets;
	        /*
			printf("[itry %d][Debug %d]route %d cur_nets.\n", itry, pid, cur_count_nets);
			printf("[itry %d][Debug %d]current net No. is: ", itry, pid);
			for (i = 0; i < cur_count_nets; i++) {
				printf("%d ", cur_nets[i]);
			}
			printf("\n");
			*/
		}

        /**********************************************Multilevel recursive routing to the crossing signals set********************************************/
		count_isignal[pid] = 0;
		for (work_id = first_work_ids[pid]; work_id < num_work_ids; work_id = 2*work_id+1) {  /* updated */
			printf("[itry %d][Debug %d]route %d work_id nets.\n", itry, pid, work_id);
		   	/* Handle the num_work_ids */	
			if (2*work_id+1 >= num_work_ids) {
				/******************Route the leaves to cur_nets****************/
				printf("[itry %d][Debug %d]route %d work_id leave nets.\n", itry, pid, work_id);
				inet = 0;
				for (k = 0; k < cur_count_nets; k++ ) {
					/* How to determine unrouted nets ? */
					inet = cur_nets[k];
					if (clb_net[inet].is_global == FALSE) {
						if(pid != 0){
							record_current_instance_each_net_tag[count_isignal[pid]++] = inet;
						}
						pathfinder_update_one_cost(trace_head[inet], -1, pres_fac);
						is_routable = breadth_first_route_net(inet, router_opts.bend_cost);
						if (!is_routable) {
							vpr_printf(TIO_MESSAGE_INFO, "Routing failed.\n");
							return (FALSE);
						}
						pathfinder_update_one_cost(trace_head[inet], 1, pres_fac);
					}
				}
				/* Finish the overall routing in own subset in own instance. */
				printf("[itry %d][Debug %d]NO.%d instance's leave set is breaking...\n", itry, pid, pid);
				break;
			}
	        assert(2*work_id+2 < num_work_ids);	/*Guojie: make it safe */
			printf("[itry %d][Debug %d]partition %d median position and %d coordianation.\n", itry, pid, median_cutline[work_id], dimension[work_id]);
			
			multilevel_paritioning_nets(cur_nets, cur_count_nets, left, &count_left, cross, &count_cross, right, &count_right, median_cutline[work_id], dimension[work_id]); 
		
			/*
			printf("[itry %d][Debug %d] %d count_left  No. net: ", itry, pid, count_left);
			for (i = 0; i < count_left; i++) {
				printf("%d ", left[i]);
			}
            printf("\n");
			printf("[itry %d][Debug %d] %d count_cross No. net: ", itry, pid, count_cross);
			for (i = 0; i < count_cross; i++) {
				printf("%d ", cross[i]);
			}
			printf("\n");
			printf("[itry %d][Debug %d] %d count_right No. net: ", itry, pid, count_right);
			for (i = 0; i < count_right; i++) {
				printf("%d ", right[i]);
			}
			printf("\n");
			*/
			//delta cost is oscillation: check_extern_nets influenced number of nets about partitioning method.
            printf("[itry %d][Debud %d]finishing %d work_id partition...\n", itry, pid, work_id);
			//count_isignal[pid] = 0;

			/* resegment the cross set which own count_cross to yield subleft, subcross, subright; */
            /* we refer the previous partitioning direction to determine this segmentation again */

		    count_subleft  = 0;
		    count_subright = 0;
		    count_subcross = 0;	
			multilevel_paritioning_nets(cross, count_cross, subleft, &count_subleft, subcross, &count_subcross, subright, &count_subright, submedian_cutline[work_id], subdimension[work_id]);
            /*
			printf("[itry %d][Debug %d] %d count_subleft  No. net: ", itry, pid, count_subleft);
			for (i = 0; i < count_subleft; i++) {
				printf("%d ", subleft[i]);
			}
            printf("\n");
			printf("[itry %d][Debug %d] %d count_subcross No. net: ", itry, pid, count_subcross);
			for (i = 0; i < count_subcross; i++) {
				printf("%d ", subcross[i]);
			}
			printf("\n");
			printf("[itry %d][Debug %d] %d count_subright No. net: ", itry, pid, count_subright);
			for (i = 0; i < count_subright; i++) {
				printf("%d ", subright[i]);
			}
			printf("\n");
			*/
			/* we plan to combine the subcrossing and subleft in same subarray */
            count_subcross_subleft_array = 0;
			for (i = 0; i < count_subcross; i++) {
				subcross_subleft_array[count_subcross_subleft_array++] = subcross[i];
			}
            for (j = 0; j < count_subleft; j++) {
				subcross_subleft_array[count_subcross_subleft_array++] = subleft[j];
			}			
            /*
			printf("[itry %d][Debug %d] %d count_subcross_subleft_array No. net: ", itry, pid, count_subcross_subleft_array);
			for (i = 0; i < count_subcross_subleft_array; i++) {
				printf("%d ", subcross_subleft_array[i]);
			}
			printf("\n");
			*/
			/************************************** Route all of subcrossing nets and subleft nets in same strategy **************************************************/
			inet = 0;
			count_routing_inet_inode = 0;
			for (k = 0; k < count_subcross_subleft_array; k++ ) {
				inet = subcross_subleft_array[k];
	            /** Handle those nets using same strategy **/			
				if (clb_net[inet].is_global == FALSE) {	

		            /* We should look this array to ensure that inet is errorless */			
					if(pid != 0){
						record_current_instance_each_net_tag[count_isignal[pid]++] = inet;
					}

					pathfinder_update_one_cost(trace_head[inet], -1, pres_fac);
					record_routing_inet_inode_item(trace_head[inet], record_before_routing_inet_inode);
					count_before_routing_inet_inode = predict_used_inode(trace_head[inet]);

					is_routable = breadth_first_route_net(inet, router_opts.bend_cost);
					/* Impossible to route? (disconnected rr_graph) */
					if (!is_routable) {
						//printf("[Debug %d]# %d net in %d work id. \n", pid, inet, work_id);	
						vpr_printf(TIO_MESSAGE_INFO, "Routing failed.\n");
						return (FALSE);
					}

					pathfinder_update_one_cost(trace_head[inet], 1, pres_fac);
					record_routing_inet_inode_item(trace_head[inet], record_after_routing_inet_inode);
					count_after_routing_inet_inode = predict_used_inode(trace_head[inet]);
				}

                qsort(record_before_routing_inet_inode, count_before_routing_inet_inode, sizeof(int), comparefunc);
                qsort(record_after_routing_inet_inode, count_after_routing_inet_inode, sizeof(int), comparefunc);
			    
				/*    
				printf("[itry %d][Debug %d] %d net node item before routing: ", itry, pid, inet);
				for (i = 0; i < count_before_routing_inet_inode; i++) {
					printf("%d ", record_before_routing_inet_inode[i]);
				}
				printf("\n");
				printf("[itry %d][Debug %d] %d net node item after routing: ", itry, pid, inet);
				for (j = 0; j < count_after_routing_inet_inode; j++) {
					printf("%d ", record_after_routing_inet_inode[j]);
				}
				printf("\n");
                */

			    before_node = 0;
			    after_node  = 0;
			    count_interleaved_routing_inet_inode = 0;
				while (1) {

					/* updated */
					if (before_node == count_before_routing_inet_inode) {
						for (; after_node < count_after_routing_inet_inode; after_node++) {
							record_interleaved_routing_inet_inode[count_interleaved_routing_inet_inode] = record_after_routing_inet_inode[after_node];
							record_interleaved_delta_inode_occ[count_interleaved_routing_inet_inode++]  = 1;				
						}
						break;
					}
					if (after_node == count_after_routing_inet_inode) {
						for (; before_node < count_before_routing_inet_inode; before_node++) {
							record_interleaved_routing_inet_inode[count_interleaved_routing_inet_inode] = record_before_routing_inet_inode[before_node];
							record_interleaved_delta_inode_occ[count_interleaved_routing_inet_inode++]  = -1;
						}
						break;
					}

					if (record_before_routing_inet_inode[before_node] < record_after_routing_inet_inode[after_node]) {
						record_interleaved_routing_inet_inode[count_interleaved_routing_inet_inode] = record_before_routing_inet_inode[before_node];
						record_interleaved_delta_inode_occ[count_interleaved_routing_inet_inode++]  = -1;
						before_node = before_node + 1;
					}
					else if (record_before_routing_inet_inode[before_node] == record_after_routing_inet_inode[after_node]) {
						before_node = before_node + 1;
						after_node  = after_node  + 1;
					}
					else { /* (record_before_routing_inet_inode[before_node] > record_after_routing_inet_inode[after_node]) */
						record_interleaved_routing_inet_inode[count_interleaved_routing_inet_inode] = record_after_routing_inet_inode[after_node];
						record_interleaved_delta_inode_occ[count_interleaved_routing_inet_inode++]  = 1;
						after_node = after_node + 1;
					}

				}
				/*
				   for (k = 0; k < count_interleaved_routing_inet_inode; k++) {
				   printf("%d[%d] ", record_interleaved_routing_inet_inode[k], record_interleaved_delta_inode_occ[k]);
				   }
				   printf("\n");
				   */
				for (interleaved_node = 0; interleaved_node < count_interleaved_routing_inet_inode; interleaved_node++) {
					store_routing_inet_inode[count_routing_inet_inode] = record_interleaved_routing_inet_inode[interleaved_node];
					store_delta_inet_inode[count_routing_inet_inode++] = record_interleaved_delta_inode_occ[interleaved_node];
				}
			}
			printf("[itry %d][Debug %d]achieving cross nets routing...\n", itry, pid);
            /*			
			printf("[itry %d][Debug %d]count_routing_inet_inode = %d \n", itry, pid, count_routing_inet_inode);
			for (i = 0; i < count_routing_inet_inode; i++) {
				printf("%d[%d] ", store_routing_inet_inode[i], store_delta_inet_inode[i]);
			}
			printf("\n");
			*/
            assert(count_routing_inet_inode < num_rr_nodes);
			/************************************************* End to route crossing nets ***************************************************************************/
			/* Determines the number of instances to the level of recursive */
			right_pid = owner_pid[2*work_id+2];
   
			printf("[itry %d][Debug %d]start to send right cost to %dth intance.\n", itry, pid, right_pid);
			MPI_Send(&count_routing_inet_inode, 1, MPI_INT, right_pid, 50, MPI_COMM_WORLD);
			MPI_Send(store_routing_inet_inode, count_routing_inet_inode, MPI_INT, right_pid, 51, MPI_COMM_WORLD);
            MPI_Send(store_delta_inet_inode, count_routing_inet_inode, MPI_INT, right_pid, 52, MPI_COMM_WORLD);
			
			/*
			for (i = 0; i < count_routing_inet_inode; i++) {
				printf("%d[%d] ", store_routing_inet_inode[i], store_delta_inet_inode[i]);
			}
			printf("\n");
            */

            /***we plan to combine the subright and right nets in same array***/
			count_subright_right_array = 0;
			for (i = 0; i < count_subright; i++) {
				subright_right_array[count_subright_right_array++] = subright[i];
			}
			for (j = 0; j < count_right; j++) {
				subright_right_array[count_subright_right_array++] = right[j];
			}

			MPI_Send(&count_subright_right_array, 1, MPI_INT, right_pid, 100, MPI_COMM_WORLD);
			MPI_Send(subright_right_array, count_subright_right_array, MPI_INT, right_pid, 101, MPI_COMM_WORLD);
			printf("[itry %d][Debug %d]send right nets to %dth intance.\n", itry, pid, right_pid);
			
			/* here is used to update the next nets information */
			cur_nets = left;
			cur_count_nets = count_left;
		    /*
			printf("[itry %d][Debug %d]cur nets = left net No. is: ", itry, pid);
			for (i = 0; i < cur_count_nets; i++) {
				printf("%d ", cur_nets[i]);
			}
			printf("\n");
			*/
		}
        printf("[itry %d][Debug %d]------------Finish Routing------------\n", itry, pid);
       
        /****************************************************************************************/
        /****************************** Parallel router completed overall signals in multilevel recursive partitioning method ***************************************/
        /*********************************************** Start to combine other instance with master instance *******************************************************/
		if (pid == 0) {	
			/** Receive used to combine with master instance such that quit the program **/
		    printf("[itry %d][Debug %d]start to wait for receiving...\n", itry, pid);
			int another_pid = 1;
			for (another_pid = 1; another_pid < num_pids; another_pid++) {

				MPI_Recv(&current_instance_total_count_nets, 1, MPI_INT, another_pid, 150, MPI_COMM_WORLD, &status_150); /* Guojie: the same issue as line above */
			    
				recv_current_instance_total_count_nets = current_instance_total_count_nets;				
				MPI_Recv(record_current_instance_each_net_tag, recv_current_instance_total_count_nets, MPI_INT, another_pid, 151, MPI_COMM_WORLD, &status_151);
				MPI_Recv(record_current_instance_count_net_nodes, recv_current_instance_total_count_nets, MPI_INT, another_pid, 152, MPI_COMM_WORLD, &status_152);
            
				/*
				printf("[itry %d][Debug %d]receive instance [%d] NO. net: ", itry, pid, another_pid);
			    for (i = 0; i < recv_current_instance_total_count_nets; i++) {
					printf("%d ", record_current_instance_each_net_tag[i]);
				}
			    printf("\n");

				printf("[itry %d][Debug %d]receive instance [%d] nodes: ", itry, pid, another_pid);
			    for (i = 0; i < recv_current_instance_total_count_nets; i++) {
					printf("%d ", record_current_instance_count_net_nodes[i]);
				}
			    printf("\n");
                */

				recv_current_instance_sum_nodes = 0;
				for (recv_current_instance_net = 0; recv_current_instance_net < recv_current_instance_total_count_nets; recv_current_instance_net++) {
					recv_current_instance_sum_nodes += record_current_instance_count_net_nodes[recv_current_instance_net];
				}
				printf("[itry %d][Debug %d]receive instance [%d] sum nodes are %d.\n", itry, pid, another_pid, recv_current_instance_sum_nodes);
				MPI_Recv(new_transfer, sizeof(t_route)*recv_current_instance_sum_nodes, MPI_BYTE, another_pid, 153, MPI_COMM_WORLD, &status_153);

				recv_current_instance_net = 0;
				recv_current_instance_total_count_nodes = 0;                    /* updated */
				for (recv_current_instance_net = 0; recv_current_instance_net < recv_current_instance_total_count_nets; recv_current_instance_net++) {
					
					inet = record_current_instance_each_net_tag[recv_current_instance_net]; /* updated */

					if (clb_net[inet].is_global == FALSE) {  /* Guojie: Free all existing trace_head[]->next and so on, otherwise you have memory that is allocated but never free */
						trace_head[inet] = (s_trace *)malloc(sizeof(s_trace)); /* Guojie: potential memory leaks--How to avoid it?  */
						recv_new_tptr = trace_head[inet];
						for (current_instance_count_net_nodes = 0; current_instance_count_net_nodes < record_current_instance_count_net_nodes[recv_current_instance_net]; current_instance_count_net_nodes++){
							recv_new_tptr->index        = new_transfer[recv_current_instance_total_count_nodes].route_index;
							recv_new_tptr->iswitch      = new_transfer[recv_current_instance_total_count_nodes].route_iswitch;
							recv_new_tptr->iblock       = new_transfer[recv_current_instance_total_count_nodes].route_iblock;
							recv_new_tptr->num_siblings = new_transfer[recv_current_instance_total_count_nodes].route_num_siblings;

							if (current_instance_count_net_nodes != record_current_instance_count_net_nodes[recv_current_instance_net]-1) {
								recv_new_tptr->next = (s_trace *)malloc(sizeof(s_trace));
							}
							else {
								recv_new_tptr->next = NULL;
							}

							recv_new_tptr = recv_new_tptr->next;
							recv_current_instance_total_count_nodes++;
						}
					}

				}
				printf("[itry %d][Debug %d]receive instance [%d] total count nodes are %d. \n", itry, pid, another_pid, recv_current_instance_total_count_nodes);
				assert(recv_current_instance_total_count_nodes == recv_current_instance_sum_nodes);
			}

			/* Make sure any CLB OPINs used up by subblocks being hooked directlyto them are reserved for that purpose. */
			if (itry == 1)
				rip_up_local_opins = FALSE;
			else
				rip_up_local_opins = TRUE;
			reserve_locally_used_opins(pres_fac, rip_up_local_opins, clb_opins_used_locally);	
			/* Attempt to certificate the feasible routing */
			success = feasible_routing();
			if (success) {
				is_success = 1;
				/* updated */
				printf("[itry %d][Debug %d]the router is successfully quiting after %d iterations...\n", itry, pid, itry);
				//return (TRUE);
			}

			printf("[itry %d][Debug %d]is_success value is %d.\n", itry, pid, is_success);
			if (is_success == 0) {
				/* Update the present cost */
				if (itry == 1)
					pres_fac = router_opts.initial_pres_fac;
				else
					pres_fac *= router_opts.pres_fac_mult;
				pres_fac = std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));
				/* Update the history and present cost */
				pathfinder_update_cost(pres_fac, router_opts.acc_fac);
				printf("[itry %d][Debug %d]to be continue iteration...\n", itry, pid);
			}
		}
		else {
     		/****************************** Reference the record_net_instance_id[inet] = pid to Collect each signal information *****************************************/
			/** Send relevant routing result **/
			current_instance_sum_nodes = 0;
			current_instance_total_count_nodes = 0;
			current_instance_total_count_nets = 0;

			current_instance_total_num_signals = count_isignal[pid]; /* updated */
	        
			/*
			printf("[itry %d][Debug %d]collection send NO. net is: ", itry, pid);
            for (i = 0; i < current_instance_total_num_signals; i++){
				printf("%d ", record_current_instance_each_net_tag[i]);
			}
            printf("\n");
            */

			for (isignal = 0; isignal < current_instance_total_num_signals; isignal++) {
				inet = record_current_instance_each_net_tag[isignal];
				current_instance_count_net_nodes = 0;
				if (clb_net[inet].is_global == FALSE) {
			    //	printf("[itry %d][Debug %d]execute copy No. %d net each node.\n", itry, pid, inet);
					new_tptr = trace_head[inet];
					for (; ;) {
						new_transfer[current_instance_total_count_nodes].route_index        = new_tptr->index;
						new_transfer[current_instance_total_count_nodes].route_iswitch      = new_tptr->iswitch;
						new_transfer[current_instance_total_count_nodes].route_iblock       = new_tptr->iblock;
						new_transfer[current_instance_total_count_nodes].route_num_siblings = new_tptr->num_siblings;

						current_instance_count_net_nodes++;
						current_instance_total_count_nodes++;

						new_tptr = new_tptr->next;
						if (new_tptr == NULL)
							break;
					}
					record_current_instance_count_net_nodes[isignal] = current_instance_count_net_nodes;
					current_instance_total_count_nets++;
				}
			}
            assert(current_instance_total_count_nets == current_instance_total_num_signals);
            
			MPI_Send(&current_instance_total_count_nets, 1, MPI_INT, 0, 150, MPI_COMM_WORLD);
			MPI_Send(record_current_instance_each_net_tag, current_instance_total_count_nets, MPI_INT, 0, 151, MPI_COMM_WORLD);
			MPI_Send(record_current_instance_count_net_nodes, current_instance_total_count_nets, MPI_INT, 0, 152, MPI_COMM_WORLD);

			printf("[itry %d][Debug %d]current_instance_total_count_nodes = %d \n", itry, pid, current_instance_total_count_nodes);
			for (current_instance_net = 0; current_instance_net < current_instance_total_count_nets; current_instance_net++) {
				current_instance_sum_nodes += record_current_instance_count_net_nodes[current_instance_net];
			}
			printf("[itry %d][Debug %d]current_instance_sum_nodes = %d \n", itry, pid, current_instance_sum_nodes);
			assert(current_instance_sum_nodes == current_instance_total_count_nodes);
			MPI_Send(new_transfer, sizeof(t_route)*current_instance_sum_nodes, MPI_BYTE, 0, 153, MPI_COMM_WORLD);

		}
	}

    free(new_transfer);

	free(overall_signals);
	free(left);
	free(right);
	free(cross);

	free(subleft);
	free(subright);
	free(subcross);
	free(subcross_subleft_array);
	free(subright_right_array);
	
	free(record_before_routing_inet_inode);      
	free(record_after_routing_inet_inode);       
	free(record_interleaved_routing_inet_inode); 
	free(record_interleaved_delta_inode_occ); 
	
	free(store_routing_inet_inode);
	free(store_delta_inet_inode);

    free(count_isignal);
    free(record_current_instance_each_net_tag);

    free(updated_overall_routing_resource_graph_cost);

	vpr_printf(TIO_MESSAGE_INFO, "Routing failed.\n");
	return (FALSE);
}

static boolean breadth_first_route_net(int inet, float bend_cost) {

	/* Uses a maze routing (Dijkstra's) algorithm to route a net.  The net       *
	 * begins at the net output, and expands outward until it hits a target      *
	 * pin.  The algorithm is then restarted with the entire first wire segment  *
	 * included as part of the source this time.  For an n-pin net, the maze     *
	 * router is invoked n-1 times to complete all the connections.  Inet is     *
	 * the index of the net to be routed.  Bends are penalized by bend_cost      *
	 * (which is typically zero for detailed routing and nonzero only for global *
	 * routing), since global routes with lots of bends are tougher to detailed  *
	 * route (using a detailed router like SEGA).                                *
	 * If this routine finds that a net *cannot* be connected (due to a complete *
	 * lack of potential paths, rather than congestion), it returns FALSE, as    *
	 * routing is impossible on this architecture.  Otherwise it returns TRUE.   */

	int i, inode, prev_node, remaining_connections_to_sink;
	float pcost, new_pcost;
	struct s_heap *current;
	struct s_trace *tptr;

	free_traceback(inet);
	breadth_first_add_source_to_heap(inet);
	mark_ends(inet);

	tptr = NULL;
	remaining_connections_to_sink = 0;

	for (i = 1; i <= clb_net[inet].num_sinks; i++) { /* Need n-1 wires to connect n pins */
		breadth_first_expand_trace_segment(tptr, remaining_connections_to_sink);
		current = get_heap_head();

		if (current == NULL) { /* Infeasible routing.  No possible path for net. */
			vpr_printf (TIO_MESSAGE_INFO, "Cannot route net #%d (%s) to sink #%d -- no possible path.\n",
					inet, clb_net[inet].name, i);
			reset_path_costs(); /* Clean up before leaving. */
			return (FALSE);
		}

		inode = current->index;

		while (rr_node_route_inf[inode].target_flag == 0) {
			pcost = rr_node_route_inf[inode].path_cost;
			new_pcost = current->cost;
			if (pcost > new_pcost) { /* New path is lowest cost. */
				rr_node_route_inf[inode].path_cost = new_pcost;
				prev_node = current->u.prev_node;
				rr_node_route_inf[inode].prev_node = prev_node;
				rr_node_route_inf[inode].prev_edge = current->prev_edge;

				if (pcost > 0.99 * HUGE_POSITIVE_FLOAT) /* First time touched. */
					add_to_mod_list(&rr_node_route_inf[inode].path_cost);

				breadth_first_expand_neighbours(inode, new_pcost, inet,
						bend_cost);
			}

			free_heap_data(current);
			current = get_heap_head();

			if (current == NULL) { /* Impossible routing. No path for net. */
				vpr_printf (TIO_MESSAGE_INFO, "Cannot route net #%d (%s) to sink #%d -- no possible path.\n",
						inet, clb_net[inet].name, i);
				reset_path_costs();
				return (FALSE);
			}

			inode = current->index;
		}

		rr_node_route_inf[inode].target_flag--; /* Connected to this SINK. */
		remaining_connections_to_sink = rr_node_route_inf[inode].target_flag;
		tptr = update_traceback(current, inet);
		free_heap_data(current);
	}

	empty_heap();
	reset_path_costs();
	return (TRUE);
}

static void breadth_first_expand_trace_segment(struct s_trace *start_ptr,
		int remaining_connections_to_sink) {

	/* Adds all the rr_nodes in the traceback segment starting at tptr (and     *
	 * continuing to the end of the traceback) to the heap with a cost of zero. *
	 * This allows expansion to begin from the existing wiring.  The            *
	 * remaining_connections_to_sink value is 0 if the route segment ending     *
	 * at this location is the last one to connect to the SINK ending the route *
	 * segment.  This is the usual case.  If it is not the last connection this *
	 * net must make to this SINK, I have a hack to ensure the next connection  *
	 * to this SINK goes through a different IPIN.  Without this hack, the      *
	 * router would always put all the connections from this net to this SINK   *
	 * through the same IPIN.  With LUTs or cluster-based logic blocks, you     *
	 * should never have a net connecting to two logically-equivalent pins on   *
	 * the same logic block, so the hack will never execute.  If your logic     *
	 * block is an and-gate, however, nets might connect to two and-inputs on   *
	 * the same logic block, and since the and-inputs are logically-equivalent, *
	 * this means two connections to the same SINK.                             */

	struct s_trace *tptr, *next_ptr;
	int inode, sink_node, last_ipin_node;

	tptr = start_ptr;
	if(tptr != NULL && rr_node[tptr->index].type == SINK) {
		/* During logical equivalence case, only use one opin */
		tptr = tptr->next;
	}

	if (remaining_connections_to_sink == 0) { /* Usual case. */
		while (tptr != NULL) {
			node_to_heap(tptr->index, 0., NO_PREVIOUS, NO_PREVIOUS, OPEN, OPEN);
			tptr = tptr->next;
		}
	}

	else { /* This case never executes for most logic blocks. */

		/* Weird case.  Lots of hacks. The cleanest way to do this would be to empty *
		 * the heap, update the congestion due to the partially-completed route, put *
		 * the whole route so far (excluding IPINs and SINKs) on the heap with cost  *
		 * 0., and expand till you hit the next SINK.  That would be slow, so I      *
		 * do some hacks to enable incremental wavefront expansion instead.          */

		if (tptr == NULL)
			return; /* No route yet */

		next_ptr = tptr->next;
		last_ipin_node = OPEN; /* Stops compiler from complaining. */

		/* Can't put last SINK on heap with NO_PREVIOUS, etc, since that won't let  *
		 * us reach it again.  Instead, leave the last traceback element (SINK) off *
		 * the heap.                                                                */

		while (next_ptr != NULL) {
			inode = tptr->index;
			node_to_heap(inode, 0., NO_PREVIOUS, NO_PREVIOUS, OPEN, OPEN);

			if (rr_node[inode].type == IPIN)
				last_ipin_node = inode;

			tptr = next_ptr;
			next_ptr = tptr->next;
		}

		/* This will stop the IPIN node used to get to this SINK from being         *
		 * reexpanded for the remainder of this net's routing.  This will make us   *
		 * hook up more IPINs to this SINK (which is what we want).  If IPIN        *
		 * doglegs are allowed in the graph, we won't be able to use this IPIN to   *
		 * do a dogleg, since it won't be re-expanded.  Shouldn't be a big problem. */

		rr_node_route_inf[last_ipin_node].path_cost = -HUGE_POSITIVE_FLOAT;

		/* Also need to mark the SINK as having high cost, so another connection can *
		 * be made to it.                                                            */

		sink_node = tptr->index;
		rr_node_route_inf[sink_node].path_cost = HUGE_POSITIVE_FLOAT;

		/* Finally, I need to remove any pending connections to this SINK via the    *
		 * IPIN I just used (since they would result in congestion).  Scan through   *
		 * the heap to do this.                                                      */

		invalidate_heap_entries(sink_node, last_ipin_node);
	}
}

static void breadth_first_expand_neighbours(int inode, float pcost, int inet,
		float bend_cost) {

	/* Puts all the rr_nodes adjacent to inode on the heap.  rr_nodes outside   *
	 * the expanded bounding box specified in route_bb are not added to the     *
	 * heap.  pcost is the path_cost to get to inode.                           */

	int iconn, to_node, num_edges;
	t_rr_type from_type, to_type;
	float tot_cost;

	num_edges = rr_node[inode].num_edges;
	for (iconn = 0; iconn < num_edges; iconn++) {
		to_node = rr_node[inode].edges[iconn];

		if (rr_node[to_node].xhigh < route_bb[inet].xmin
				|| rr_node[to_node].xlow > route_bb[inet].xmax
				|| rr_node[to_node].yhigh < route_bb[inet].ymin
				|| rr_node[to_node].ylow > route_bb[inet].ymax)
			continue; /* Node is outside (expanded) bounding box. */

		tot_cost = pcost + get_rr_cong_cost(to_node);

		if (bend_cost != 0.) {
			from_type = rr_node[inode].type;
			to_type = rr_node[to_node].type;
			if ((from_type == CHANX && to_type == CHANY)
					|| (from_type == CHANY && to_type == CHANX))
				tot_cost += bend_cost;
		}

		node_to_heap(to_node, tot_cost, inode, iconn, OPEN, OPEN);
	}
}

static void breadth_first_add_source_to_heap(int inet) {

	/* Adds the SOURCE of this net to the heap.  Used to start a net's routing. */

	int inode;
	float cost;

	inode = net_rr_terminals[inet][0]; /* SOURCE */
	cost = get_rr_cong_cost(inode);

	node_to_heap(inode, cost, NO_PREVIOUS, NO_PREVIOUS, OPEN, OPEN);
}

static void multilevel_paritioning_nets(int *nets, int count_nets, int *left, int *count_left, int *cross, int *count_cross, int *right, int *count_right, int med_cutline, int dim ) {

	int pid;
   	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    //printf("--------------- This multilevel partitioning from %dth intance ---------------\n", pid);
	/* Here perform the partitioning nets. */
	int a_net = 0;
	int i = 0, j = 0, k = 0;

    int inet = 0;
    int	jpin = 0, kpin = 0;
    int crossing = 0, left_down = 0, right_up = 0;

	int *flag_nets = (int *)malloc(sizeof(int)*num_nets);
	assert(flag_nets);
	memset(flag_nets, 0, sizeof(int)*num_nets);

    /* 	
    printf("\n[Debug %d]subfunction included No. net: ", pid);	
	for (a_net = 0; a_net < count_nets; a_net++) {
		printf("%d ", nets[a_net]);
	}
	printf("\n");
	*/

	/* Make these all nets divided into left, cross and right */
	for (a_net = 0; a_net < count_nets; a_net++) {
		inet = nets[a_net];
		if (clb_net[inet].is_global == FALSE) {
			int sum_sinks = clb_net[inet].num_sinks;
			int *flag_pins = (int *)malloc(sizeof(int)*(sum_sinks+1));
			assert(flag_pins);
			memset(flag_pins, 0, sizeof(int)*(sum_sinks+1));
			/* flag to each terminials of signal to vertical direction. left = 1 and right = 2 */
		    /* 0: x-coordination  1: y-coordination */
			/* the med_cutline includes position to cutlines. */
			if (dim == 0) { 
				for (jpin = 0; jpin <= sum_sinks; jpin++) {
					if (block[clb_net[inet].node_block[jpin]].x < med_cutline) {
						flag_pins[jpin] = 1; // left = 1
					}
					else {
						flag_pins[jpin] = 2; // right = 2
					}
				}
			}
			else {
				for (jpin = 0; jpin <= sum_sinks; jpin++) {
					if (block[clb_net[inet].node_block[jpin]].y < med_cutline) {
						flag_pins[jpin] = 1; // down = 1
					}
					else {
						flag_pins[jpin] = 2; // up   = 2
					}
				}
			}
			/* Starts to classify these nets, here is core partition spirits. */
			for (kpin = 1; kpin <= sum_sinks; kpin++) {
				if (flag_pins[kpin] != flag_pins[kpin-1]) {
					flag_nets[inet] = 3;
					crossing++;
					break;
				}
				else if (flag_pins[kpin] == 1) {
					flag_nets[inet] = flag_pins[kpin];
					left_down++;
					break;
				}
				else { /* (flag_pins[kpin] == 2)*/
					flag_nets[inet] = flag_pins[kpin];
					right_up++;
					break;
				}
			} 
			free(flag_pins);
		}
	}
    //	printf("left_down = %d crossing = %d right_up = %d from %dth intance \n", left_down, crossing, right_up, pid);
	/* Allocate the new array to record this data such that recursive routing. */
	i = 0, j = 0, k = 0;
	for (a_net = 0; a_net < count_nets; a_net++) {
		inet = nets[a_net];
		if (clb_net[inet].is_global == FALSE) {
			if (flag_nets[inet] == 1) { 
				left[i++]  = inet;
				assert(i <= left_down);
			}
			if (flag_nets[inet] == 2) {
				right[j++] = inet;
				assert(j <= right_up);
			}
			if (flag_nets[inet] == 3) {
				cross[k++] = inet;
				assert(k <= crossing);
			}
		}
	}

	if (flag_nets) {
		free(flag_nets);
	}

  /* 
	int ileft, jright, kcross;
	printf("[Debug %d]------------------left %d nets---------------------\n", pid, i);
	for (ileft = 0; ileft < i; ileft++) {
		printf("%d ", left[ileft]);
	}
	printf("\n");

	printf("[Debug %d]------------------right %d nets---------------------\n", pid, j);
	for (jright = 0; jright < j; jright++ ) {
		printf("%d ", right[jright]);
	}
	printf("\n");

	printf("[Debug %d]------------------cross %d nets---------------------\n", pid, k);
	for (kcross = 0; kcross < k; kcross++) {
		printf("%d ", cross[kcross]);
	}
	printf("\n");
 */
    *count_left  = i;
	*count_right = j;
	*count_cross = k;
//	printf("[Debug %d]subfunction count_left = %d count_cross = %d count_right = %d\n", pid, *count_left, *count_cross, *count_right);
	//free(flag_nets);
//	printf("[Debug %d]subfunction is completely quiting...\n", pid);
}

static void record_routing_inet_inode_item(struct s_trace *routing_net_start, int *routing_inet_inode)
{
	struct s_trace *tptr;
	int i = 0, inode = 0;
	tptr = routing_net_start;
	
	if (tptr == NULL)
		return;

	for (;;) {
		inode = tptr->index;
		routing_inet_inode[i] = inode;
		++i;
		if (rr_node[inode].type == SINK) {
			tptr = tptr->next;
			if (tptr == NULL)
				break;
		}
		tptr = tptr->next;
	}
}

static int predict_used_inode(struct s_trace *routing_segment_start)
{
	struct s_trace *tptr;
	int i = 0, inode = 0;

	tptr = routing_segment_start;
	if (tptr == NULL)
	   return 1;

    for (;;) {
		inode = tptr->index;
		++i;
		if (rr_node[inode].type == SINK) {
			tptr = tptr->next;
			if (tptr == NULL)
			   return i;
		}
		tptr = tptr->next;
	}	
}

int comparefunc (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

