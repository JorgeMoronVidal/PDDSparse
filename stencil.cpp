#include "stencil.hpp"
Stencil::Stencil(void){
    //All variables of the class are emptied
    pos_north.clear(); index_north.clear();
    pos_south.clear(); index_south.clear();
    pos_east.clear(); index_east.clear();
    pos_west.clear(); index_west.clear();
}
void Stencil::Init(int dir,
              bool interior, 
              std::vector<int> current_index,
              std::vector<Eigen::VectorXf> & node_position,
              std::vector<std::vector<int> > & node_index,
              std::map<std::vector<int>, int> & interface_map,
              std::vector<int> & number_interfaces, 
              float *g_parameters){
	std::vector< std::vector<int> > vvaux_north, vvaux_south, vvaux_east, vvaux_west;
    std::vector< int > vaux;
    vaux.resize(2);
    for(uint16_t i = 0; i < 4; i++) stencil_parameters[i] = global_parameters[i] = g_parameters[i];
    if(interior == true){
        switch (dir) {
    
            case 0: 
                //Spline is Horizontal

                //Top spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 2;
                stencil_parameters[3] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_north.push_back(vaux);
                

                //Bottom spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 2;
                stencil_parameters[1] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_south.push_back(vaux);

                //Letf-down spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] - 1;
                stencil_parameters[0] = node_position[node_index[interface_map[vaux]].back()][0];
                vvaux_west.push_back(vaux);

                //Left-up spline interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] + 1;
                vvaux_west.push_back(vaux);

                //Right-down spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 1;
                vvaux_east.push_back(vaux);

                //Right-up spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 1;
                stencil_parameters[2] = node_position[node_index[interface_map[vaux]].back()][0];
                vvaux_east.push_back(vaux);

                
            break;
    
            case 1:
                //Vertical

                //Right spline Interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1];
                stencil_parameters[2] = node_position[node_index[interface_map[vaux]].back()][0];
                vvaux_east.push_back(vaux);

                //Left spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1];
                stencil_parameters[0] = node_position[node_index[interface_map[vaux]].back()][0];
                vvaux_west.push_back(vaux);

                //Letf-up spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 1;
                vvaux_north.push_back(vaux);

                //Right-up spline interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1] + 1;
                stencil_parameters[3] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_north.push_back(vaux);

                //Letf-down spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 1;
                stencil_parameters[1] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_south.push_back(vaux);

                //Right-dowm spline interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1] - 1;
                vvaux_south.push_back(vaux);

    
            break;
    
            default :
            std::cout <<"Wrong direction parameter. Stencils won't be properly built."<< '\n';
            }

    }else{
        switch (dir) {
    
            case 0: 
                //Horizontal

                //Top spline Interface
                if(current_index[1] + 2 <= 2*number_interfaces[1] - 3){
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 2;
                stencil_parameters[3] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_north.push_back(vaux);
                }

                //Bottom spline Interface
                if(current_index[1]-2 >= 0){
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 2;
                stencil_parameters[1] = node_position[node_index[interface_map[vaux]].back()][1];
                vvaux_south.push_back(vaux);
                }

                //Left splines
                if(current_index[0] > 0){
                    //Letf-down spline Interface
                    if(current_index[1]-1 >= 0){
                        vaux[0] = current_index[0] - 1;
                        vaux[1] = current_index[1] - 1;
                        stencil_parameters[0] = node_position[node_index[interface_map[vaux]].back()][0];
                        vvaux_west.push_back(vaux);
                    }
                    //Left-up spline interface
                    if(current_index[1]+1 <= 2*number_interfaces[1] - 2){
                        vaux[0] = current_index[0] - 1;
                        vaux[1] = current_index[1] + 1;
                        vvaux_west.push_back(vaux);
                    }
                }

                //Right splines
                if(current_index[0] < number_interfaces[0] - 1){
                    //Right-down spline Interface
                    if(current_index[1]-1 >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] - 1;
                        vvaux_east.push_back(vaux);
                    }
                    //Right-up spline interface
                    if(current_index[1]+1 <= 2*number_interfaces[1] - 2){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] + 1;
                        stencil_parameters[2] = node_position[node_index[interface_map[vaux]].back()][0];
                        vvaux_east.push_back(vaux);
                    }
                }

            break;
    
            case 1:
                //Vertical

                //Right spline Interface
                if(current_index[0] + 1 <= number_interfaces[0] - 2){
                    vaux[0] = current_index[0] + 1;
                    vaux[1] = current_index[1];
                    stencil_parameters[2] = node_position[node_index[interface_map[vaux]].back()][0];
                    vvaux_east.push_back(vaux);
                }

                //Left spline Interface
                if(current_index[0] - 1 >= 0){
                    vaux[0] = current_index[0] - 1;
                    vaux[1] = current_index[1];
                    stencil_parameters[0] = node_position[node_index[interface_map[vaux]].back()][0];
                    vvaux_west.push_back(vaux);
                }

                //Top splines
                if(current_index[1] + 1 <= 2*number_interfaces[1] - 3 ){
                    //Letf-up spline Interface
                    if(current_index[0] >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] + 1;
                        vvaux_north.push_back(vaux);
                    }
                    //Right-up spline interface
                    if(current_index[0] + 1 <= number_interfaces[0] - 1){
                        vaux[0] = current_index[0] + 1;
                        vaux[1] = current_index[1] + 1;
                        stencil_parameters[3] = node_position[node_index[interface_map[vaux]].back()][1];
                        vvaux_north.push_back(vaux);
                    }
                }

                //Bottom splines
                if(current_index[1] - 1 >= 0){
                    //Letf-down spline Interface
                    if(current_index[0] >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] - 1;
                        stencil_parameters[1] = node_position[node_index[interface_map[vaux]].back()][1];
                        vvaux_south.push_back(vaux);
                    }
                    //Right-dowm spline interface
                    if(current_index[0] + 1 <= number_interfaces[0] - 1){
                        vaux[0] = current_index[0] + 1;
                        vaux[1] = current_index[1] - 1;
                        vvaux_south.push_back(vaux);
                    }
                }

    
            break;
    
            default :
            std::cout <<"Wrong direction parameter. Stencils won't be properly built."<< '\n';
            } 
    
    }


    //All variables of the class are emptied
    pos_north.clear(); index_north.clear(); G_north.clear();
    pos_south.clear(); index_south.clear(); G_south.clear();
    pos_east.clear(); index_east.clear(); G_east.clear();
    pos_west.clear(); index_west.clear();G_west.clear();
    ipsi_north.resize(0,0); ipsi_south.resize(0,0);
    ipsi_east.resize(0,0); ipsi_west.resize(0,0);
    if(vvaux_north.size() > 0){
    	for(std::vector<std::vector<int> >::iterator inter = vvaux_north.begin();
        inter != vvaux_north.end();
        inter++){
    		for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
    			node != node_index[interface_map[*inter]].end();
    			node++){
    			index_north.push_back(*node);
    		}
        }
        std::vector<int>::iterator v_north = std::unique(index_north.begin(), index_north.end());
        index_north.erase(v_north, index_north.end());
        for(unsigned int i = 0; i < index_north.size(); i++) pos_north.push_back(node_position[index_north[i]]);
        for(unsigned int i = 0; i < pos_north.size(); i++) G_north.push_back(0.0f);
    }
    
    if(vvaux_south.size() > 0){
    	for(std::vector<std::vector<int> >::iterator inter = vvaux_south.begin();
        inter != vvaux_south.end();
        inter++){
    		for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
    			node != node_index[interface_map[*inter]].end();
    			node++){
    			index_south.push_back(*node);
    		}
        }
        std::vector<int>::iterator v_south = std::unique(index_south.begin(), index_south.end());
        index_south.erase(v_south, index_south.end()); 
        for(unsigned int i = 0; i < index_south.size(); i++) pos_south.push_back(node_position[index_south[i]]);
        for(unsigned int i = 0; i < pos_south.size(); i++) G_south.push_back(0.0f);
    }

    if(vvaux_east.size() > 0){
    	for(std::vector<std::vector<int> >::iterator inter = vvaux_east.begin();
        inter != vvaux_east.end();
        inter++){
    		for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
    			node != node_index[interface_map[*inter]].end();
    			node++){
    			index_east.push_back(*node);
    		}
        }
        std::vector<int>::iterator v_east = std::unique(index_east.begin(), index_east.end());
        index_east.erase(v_east, index_east.end()); 
        for(unsigned int i = 0; i < index_east.size(); i++) pos_east.push_back(node_position[index_east[i]]);
        for(unsigned int i = 0; i < pos_east.size(); i++) G_east.push_back(0.0f);
    }

    if(vvaux_west.size() > 0){
    	for(std::vector<std::vector<int> >::iterator inter = vvaux_west.begin();
        inter != vvaux_west.end();
        inter++){
    		for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
    			node != node_index[interface_map[*inter]].end();
    			node++){
    			index_west.push_back(*node);
    		}
        }
        std::vector<int>::iterator v_west = std::unique(index_west.begin(), index_west.end());
        index_west.erase(v_west, index_west.end());
        for(unsigned int i = 0; i < index_west.size(); i++) pos_west.push_back(node_position[index_west[i]]);
        for(unsigned int i = 0; i < pos_west.size(); i++) G_west.push_back(0.0f);
    }
	counter_north = 0; counter_south = 0;
	counter_east = 0; counter_west = 0;

}

void Stencil::Init(std::vector<int> interface_1,
          std::vector<int> interface_2,
          std::vector<int> interface_3,
          std::vector<int> interface_4,
          std::vector<Eigen::VectorXf> & node_position,
          std::vector<std::vector<int> > & node_index,
          std::map<std::vector<int>, int> & interface_map,
          std::vector<int> & number_interfaces, 
          float *g_parameters){
    std::vector< std::vector<int> > vvaux_north, vvaux_south, vvaux_east, vvaux_west;
    std::vector< int > vaux, interface_east, interface_west,interface_north, interface_south;
    vaux.resize(2);
    for(uint16_t i = 0; i < 4; i++) stencil_parameters[i] = global_parameters[i] = g_parameters[i];
    int max_v = 0, max_h = 0;
    max_h = std::max(interface_1[0], interface_2[0]);
    max_h = std::max(max_h, interface_3[0]);
    max_h = std::max(max_h, interface_4[0]);
    max_v = std::max(interface_1[1], interface_2[1]);
    max_v = std::max(max_v, interface_3[1]);
    max_v = std::max(max_v, interface_4[1]);
    if(interface_1[1]%2 == 1){
        //Horizontal interface
        if(interface_1[0] == max_h) interface_east = interface_1;
        else interface_west = interface_1;
    } else {
        //Vertical interface
        if(interface_1[1] == max_v) interface_north = interface_1;
        else interface_south = interface_1;
    }
    if(interface_2[1]%2 == 1){
        //Horizontal interface
        if(interface_2[0] == max_h) interface_east = interface_2;
        else interface_west = interface_2;
    } else {
        //Vertical interface
        if(interface_2[1] == max_v) interface_north = interface_2;
        else interface_south = interface_2;
    }
    if(interface_3[1]%2 == 1){
        //Horizontal interface
        if(interface_3[0] == max_h) interface_east = interface_3;
        else interface_west = interface_3;
    } else {
        //Vertical interface
        if(interface_3[1] == max_v) interface_north = interface_3;
        else interface_south = interface_3;
    }
    if(interface_4[1]%2 == 1){
        //Horizontal interface
        if(interface_4[0] == max_h) interface_east = interface_4;
        else interface_west = interface_4;
    } else {
        //Vertical interface
        if(interface_4[1] == max_v) interface_north = interface_4;
        else interface_south = interface_4;
    }

    //Top spline Interface
    if(interface_east[1] + 2 <= 2*number_interfaces[1] - 3){
        vaux[0] = interface_west[0];
        vaux[1] = interface_west[1] + 2;
        vvaux_north.push_back(vaux);
        vaux[0] = interface_east[0];
        vaux[1] = interface_east[1] + 2;
        stencil_parameters[3] = node_position[node_index[interface_map[vaux]].back()][1];
        vvaux_north.push_back(vaux);
       
    }
    
    //Bottom spline Interface
    if(interface_east[1]-2 >= 1){
        vaux[0] = interface_west[0];
        vaux[1] = interface_west[1] - 2;
        vvaux_south.push_back(vaux);
        vaux[0] = interface_east[0];
        vaux[1] = interface_east[1] - 2;
        stencil_parameters[1] = node_position[node_index[interface_map[vaux]].back()][1];
        vvaux_south.push_back(vaux);
       
    }

    //Left splines
    if(interface_north[0] > 0){
        vaux[0] = interface_south[0] - 1;
        vaux[1] = interface_south[1];
        stencil_parameters[0] = node_position[node_index[interface_map[vaux]].back()][0];
        vvaux_west.push_back(vaux);
        vaux[0] = interface_north[0] - 1;
        vaux[1] = interface_north[1];
        vvaux_west.push_back(vaux);       
    }

    //Right splines
    if(interface_north[0] < number_interfaces[0] - 2){
        vaux[0] = interface_south[0] + 1;
        vaux[1] = interface_south[1];
        stencil_parameters[2] = node_position[node_index[interface_map[vaux]].back()][0];
        vvaux_east.push_back(vaux);
        vaux[0] = interface_north[0] + 1;
        vaux[1] = interface_north[1];
        vvaux_east.push_back(vaux);       
    }
    //All variables of the class are emptied
    pos_north.clear(); index_north.clear(); G_north.clear();
    pos_south.clear(); index_south.clear(); G_south.clear();
    pos_east.clear(); index_east.clear(); G_east.clear();
    pos_west.clear(); index_west.clear();G_west.clear();
    ipsi_north.resize(0,0); ipsi_south.resize(0,0);
    ipsi_east.resize(0,0); ipsi_west.resize(0,0);

    if(vvaux_north.size() > 0){
        for(std::vector<std::vector<int> >::iterator inter = vvaux_north.begin();
        inter != vvaux_north.end();
        inter++){
            for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
                node != node_index[interface_map[*inter]].end();
                node++){
                index_north.push_back(*node);
            }
        }
        std::vector<int>::iterator v_north = std::unique(index_north.begin(), index_north.end());
        index_north.erase(v_north, index_north.end());
        pos_north.clear(); 
        for(unsigned int i = 0; i < index_north.size(); i++) pos_north.push_back(node_position[index_north[i]]);
        G_north.clear();
        for(unsigned int i = 0; i < pos_north.size(); i++) G_north.push_back(0.0f);
    }
    
    if(vvaux_south.size() > 0){
        for(std::vector<std::vector<int> >::iterator inter = vvaux_south.begin();
        inter != vvaux_south.end();
        inter++){
            for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
                node != node_index[interface_map[*inter]].end();
                node++){
                index_south.push_back(*node);
            }
        }
        std::vector<int>::iterator v_south = std::unique(index_south.begin(), index_south.end());
        index_south.erase(v_south, index_south.end()); 
        pos_south.clear();
        for(unsigned int i = 0; i < index_south.size(); i++) pos_south.push_back(node_position[index_south[i]]);
        G_south.clear();
        for(unsigned int i = 0; i < pos_south.size(); i++) G_south.push_back(0.0f);
    }

    if(vvaux_east.size() > 0){
        for(std::vector<std::vector<int> >::iterator inter = vvaux_east.begin();
        inter != vvaux_east.end();
        inter++){
            for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
                node != node_index[interface_map[*inter]].end();
                node++){
                index_east.push_back(*node);
            }
        }
        std::vector<int>::iterator v_east = std::unique(index_east.begin(), index_east.end());
        index_east.erase(v_east, index_east.end()); 
        pos_east.clear();
        for(unsigned int i = 0; i < index_east.size(); i++) pos_east.push_back(node_position[index_east[i]]);
        G_east.clear();
        for(unsigned int i = 0; i < pos_east.size(); i++) G_east.push_back(0.0f);
    }

    if(vvaux_west.size() > 0){
        for(std::vector<std::vector<int> >::iterator inter = vvaux_west.begin();
        inter != vvaux_west.end();
        inter++){
            for(std::vector<int>::iterator node = node_index[interface_map[*inter]].begin();
                node != node_index[interface_map[*inter]].end();
                node++){
                index_west.push_back(*node);
            }
        }
        std::vector<int>::iterator v_west = std::unique(index_west.begin(), index_west.end());
        index_west.erase(v_west, index_west.end());
        pos_west.clear(); 
        for(unsigned int i = 0; i < index_west.size(); i++) pos_west.push_back(node_position[index_west[i]]);
        G_west.clear();
        for(unsigned int i = 0; i < pos_west.size(); i++) G_west.push_back(0.0f);
    }
    counter_north = 0; counter_south = 0;
    counter_east = 0; counter_west = 0;
}
void Stencil::Reset(void){
    G_north.clear();
    for(unsigned int i = 0; i < pos_north.size(); i++) G_north.push_back(0.0f);
    G_south.clear();
    for(unsigned int i = 0; i < pos_south.size(); i++) G_south.push_back(0.0f);
    G_east.clear();
    for(unsigned int i = 0; i < pos_east.size(); i++) G_east.push_back(0.0f);
    G_west.clear();
    for(unsigned int i = 0; i < pos_west.size(); i++) G_west.push_back(0.0f);
    counter_north = 0; counter_south = 0;
    counter_east = 0; counter_west = 0;
}

Eigen::MatrixXf Stencil::Compute_ipsi(std::vector<Eigen::VectorXf> & sten_position,
                  BVP bvp,
                  float c2){
    //Auxiliary vectors 
    Eigen::VectorXf vaux, sincrement;
    //Psi Matrix, its inverse and Identity are created 
    Eigen::MatrixXd Psi, iPsi, I;
    //Matrix are resized
    Psi.resize(sten_position.size(), sten_position.size());
    I.resize(sten_position.size(), sten_position.size());
    I.setIdentity();
    //And fullfilled
    for(unsigned int i = 0; i < sten_position.size(); i ++){
        for(unsigned int j = 0; j < sten_position.size(); j ++){
            Psi(i,j) = bvp.rbf.Value(sten_position[i], sten_position[j], c2);
        }
    }
    float err, cond;
    iPsi= Psi.inverse();
    Eigen::BiCGSTAB<Eigen::MatrixXd> CGS;
    CGS.compute(Psi);
    iPsi = CGS.solveWithGuess(I,iPsi);
    err = CGS.error();
    cond = Psi.norm()*iPsi.norm();
    //std::cout << "Inverse matrix computed with error "<< err <<" and condition number" << cond << std::endl;
    Psi.resize(0,0);
    I.resize(0,0);
    Eigen::MatrixXf iPsi_return;
    iPsi_return = iPsi.cast<float>();
    return iPsi_return;
}

bool Stencil::AreSame(float a, float b)
{
    return fabs(a - b) < EPSILON;
}

void Stencil::Compute_ipsi(BVP bvp, float c2){
	if(pos_north.size() > 0){
		ipsi_north = Compute_ipsi(pos_north, bvp, c2);
	}
	if(pos_south.size() > 0){
		ipsi_south = Compute_ipsi(pos_south, bvp, c2);
	}
	if(pos_east.size() > 0){
		ipsi_east = Compute_ipsi(pos_east, bvp, c2);
	}
	if(pos_west.size() > 0){
		ipsi_west = Compute_ipsi(pos_west, bvp, c2);
	}
}

int Stencil::G_update(Eigen::VectorXf X, float Y, BVP bvp, float c2){
	float H_ij ;
	if(AreSame(X(1),stencil_parameters[1])){
		//South stencil
        //std::cout << "\t south \n";
		counter_south ++;
		for(unsigned int i = 0; i < G_south.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_south.size(); j++){
				H_ij += ipsi_south(j,i)*bvp.rbf.Value(pos_south[j],X,c2);
			}
			G_south[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(1),stencil_parameters[3])){
		//North stencil
        //std::cout << "\t north \n";
		counter_north ++;
		for(unsigned int i = 0; i < G_north.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_north.size(); j++){
				H_ij += ipsi_north(j,i)*bvp.rbf.Value(pos_north[j],X,c2);
			}
			G_north[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(0),stencil_parameters[0])){
		//West stencil
        //std::cout << "\t west \n";
		counter_west ++;
		for(unsigned int i = 0; i < G_west.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_west.size(); j++){
				H_ij += ipsi_west(j,i)*bvp.rbf.Value(pos_west[j],X,c2);
			}
			G_west[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(0),stencil_parameters[2])){
		//East stencil
        //std::cout << "\t east \n";
		counter_east ++;
		for(unsigned int i = 0; i < G_east.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_east.size(); j++){
				H_ij += ipsi_east(j,i)*bvp.rbf.Value(pos_east[j],X,c2);
			}
			G_east[i] += -H_ij*Y;
		}
		return 1;
	}
	std::cout << "Something went wrong computing H.";
	return 0;
}

float Stencil::Test_Interpolator(Eigen::VectorXf X, BVP bvp, float c2){
    float uH_ij;
    if(AreSame(X(1),stencil_parameters[1])){
        //South stencil
        //std::cout << "\t south \n";
        counter_south ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_south.size(); i++){
            for(unsigned int j = 0; j < G_south.size(); j++){
                uH_ij += ipsi_south(j,i)*bvp.rbf.Value(pos_south[j],X,c2)*bvp.u.Value(pos_south[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(1),stencil_parameters[3])){
        //North stencil
        //std::cout << "\t north \n";
        counter_north ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_north.size(); i++){
            for(unsigned int j = 0; j < G_north.size(); j++){
                uH_ij += ipsi_north(j,i)*bvp.rbf.Value(pos_north[j],X,c2)*bvp.u.Value(pos_north[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(0),stencil_parameters[0])){
        //West stencil
        //std::cout << "\t west \n";
        counter_west ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_west.size(); i++){
            for(unsigned int j = 0; j < G_west.size(); j++){
                uH_ij += ipsi_west(j,i)*bvp.rbf.Value(pos_west[j],X,c2)*bvp.u.Value(pos_west[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(0),stencil_parameters[2])){
        //East stencil
        //std::cout << "\t east \n";
        counter_east ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_east.size(); i++){
            for(unsigned int j = 0; j < G_east.size(); j++){
                uH_ij += ipsi_east(j,i)*bvp.rbf.Value(pos_east[j],X,c2)*bvp.u.Value(pos_east[i],0.0f);
            }
        }
        return uH_ij;
    }
    std::cout << "Something went wrong computing H.";
    return 0.0f;
}

void Stencil::G_return(std::vector<int> & G_j, std::vector<float> & G){
	G_j.clear();
	G.clear();
	float norm = 1.0f/(counter_north + counter_south + counter_east + counter_west);
	if(pos_north.size() > 0){
        //norm = 1.0f/(counter_north);
		for(unsigned int i = 0; i < index_north.size(); i++){
			G_j.push_back(index_north[i]);
			G.push_back(G_north[i]*norm);
		}
	}
	if(pos_south.size() > 0){
        //norm = 1.0f/(counter_south);
		for(unsigned int i = 0; i < index_south.size(); i++){
			G_j.push_back(index_south[i]);
			G.push_back(G_south[i]*norm);
		}
	}
	if(pos_east.size() > 0){
        //norm = 1.0f/(counter_east);
		for(unsigned int i = 0; i < index_east.size(); i++){
			G_j.push_back(index_east[i]);
			G.push_back(G_east[i]*norm);
		}
	}
	if(pos_west.size() > 0){
        //norm = 1.0f/(counter_west);
		for(unsigned int i = 0; i < index_west.size(); i++){
			G_j.push_back(index_west[i]);
			G.push_back(G_west[i]*norm);
		}
	}

}
void Stencil::G_return_withrep(std::vector<int> & G_j, std::vector<float> & G, int N_tray){
    G_j.clear();
    G.clear();
    float norm = 1.0f/(N_tray*1.0f);
    std::map<int,float> G_map;
    if(pos_north.size() > 0){
        for(unsigned int i = 0; i < index_north.size(); i++){
            G_map[index_north[i]] = 0.0f;
        }
    }
    if(pos_south.size() > 0){
        for(unsigned int i = 0; i < index_south.size(); i++){
            G_map[index_south[i]] = 0.0f;
        }
    }
    if(pos_east.size() > 0){
        for(unsigned int i = 0; i < index_east.size(); i++){
            G_map[index_east[i]] = 0.0f;
        }
    }
    if(pos_west.size() > 0){
        for(unsigned int i = 0; i < index_west.size(); i++){
            G_map[index_west[i]] = 0.0f;
        }
    }
    if(pos_north.size() > 0){
        for(unsigned int i = 0; i < index_north.size(); i++){
            G_map[index_north[i]] =+ G_north[i]*norm;
        }
    }
    if(pos_south.size() > 0){
        for(unsigned int i = 0; i < index_south.size(); i++){
            G_map[index_south[i]] =+ G_south[i]*norm;
        }
    }
    if(pos_east.size() > 0){
        for(unsigned int i = 0; i < index_east.size(); i++){
            G_map[index_east[i]] =+ G_east[i]*norm;
        }
    }
    if(pos_west.size() > 0){
        for(unsigned int i = 0; i < index_west.size(); i++){
            G_map[index_west[i]] =+ G_west[i]*norm;
        }
    }
    for(std::map<int, float>::iterator it = G_map.begin();
        it != G_map.end();
        it ++){
        G_j.push_back(it->first);
        G.push_back(it->second);
    }

}
int Stencil::G_Test_update(Eigen::VectorXf X){
    if(AreSame(X(1),stencil_parameters[1])){
        //South stencil
        //std::cout << "\t south \n";
        for(unsigned int i = 0; i < G_south.size(); i++) G_south[i] = (float)index_south[i];
            return 1;
    }
    if(AreSame(X(1),stencil_parameters[3])){
        //North stencil
        for(unsigned int i = 0; i < G_north.size(); i++) G_north[i] = (float) index_north[i];
            return 1;
    }
    if(AreSame(X(0),stencil_parameters[0])){
        //West stencil
        //std::cout << "\t west \n";
        for(unsigned int i = 0; i < G_west.size(); i++) G_west[i] = (float) index_west[i];
            return 1;
    }
    if(AreSame(X(0),stencil_parameters[2])){
        //East stencil
        for(unsigned int i = 0; i < G_east.size(); i++) G_east[i] = (float) index_east[i];
            return 1;
    }
    std::cout << "Something went wrong computing H.";
    return 0;
}

void Stencil::G_Test_return(std::vector<int> & stencil_index, std::vector<float> & G){
    stencil_index.clear();
    G.clear();
    if(pos_north.size() > 0){
        //norm = 1.0f/(counter_north);
        for(unsigned int i = 0; i < index_north.size(); i++){
            stencil_index.push_back(index_north[i]);
            G.push_back(G_north[i]);
        }
    }
    if(pos_south.size() > 0){
        //norm = 1.0f/(counter_south);
        for(unsigned int i = 0; i < index_south.size(); i++){
            stencil_index.push_back(index_south[i]);
            G.push_back(G_south[i]);
        }
    }
    if(pos_east.size() > 0){
        //norm = 1.0f/(counter_east);
        for(unsigned int i = 0; i < index_east.size(); i++){
            stencil_index.push_back(index_east[i]);
            G.push_back(G_east[i]);
        }
    }
    if(pos_west.size() > 0){
        //norm = 1.0f/(counter_west);
        for(unsigned int i = 0; i < index_west.size(); i++){
            stencil_index.push_back(index_west[i]);
            G.push_back(G_west[i]);
        }
    }

}

void Stencil::Print(int node_index){
    char fname[100];
    FILE *pf;
    sprintf(fname,"Output/Debug/stencil_%d.txt", node_index);
    pf = fopen(fname, "w");
    fprintf(pf,"index,x,y,sten\n");
    for(unsigned int i = 0; i < index_south.size(); i++) fprintf(pf,"%d,%.3f,%.3f,south\n",index_south[i], pos_south[i](0), pos_south[i](1));
    for(unsigned int i = 0; i < index_east.size(); i++) fprintf(pf,"%d,%.3f,%.3f,east\n",index_east[i], pos_east[i](0), pos_east[i](1));
    for(unsigned int i = 0; i < index_north.size(); i++) fprintf(pf,"%d,%.3f,%.3f,north\n",index_north[i], pos_north[i](0), pos_north[i](1));
    for(unsigned int i = 0; i < index_west.size(); i++) fprintf(pf,"%d,%.3f,%.3f,west\n",index_west[i], pos_west[i](0), pos_west[i](1));
    fclose(pf);
    sprintf(fname,"Output/Debug/boundary_stencil_%d.txt", node_index);
    pf = fopen(fname, "w");
    fprintf(pf,"x,y\n");
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[1]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[2], stencil_parameters[1]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[2], stencil_parameters[3]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[3]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[1]);

}

