//
// Created by acsr on 4/26/21.
//

#ifndef ACSR_PLANNER_NANOWIRE_PLANNER_HPP
#define ACSR_PLANNER_NANOWIRE_PLANNER_HPP
#include <utility>
#include <numeric>
#include <random>
#include <atomic>
#include "tree_node.hpp"
#include "dynamics.hpp"

namespace acsr {
    
    /***
     * abstract class for planner
     */
    template<typename DynamicsModel,typename Map>
    struct Planner {
    public:
        constexpr static size_t nx_=DynamicsModel::nx;
        constexpr static size_t nu_=DynamicsModel::nu;

    //protected:
        using StateType = std::array<double,nx_>;
        using ControlType = std::array<double,nu_>;
        using TreeNodeType = TreeNode<double,nx_,nu_>;
        using TreeNodePtr = std::shared_ptr<TreeNodeType>;
        using NodeType = Node<double,nx_>;
        using NodePtr = std::shared_ptr<NodeType>;
        
    protected:
        std::shared_ptr<DynamicsModel> dynamics_model_ptr_; ///nanowire system
        std::shared_ptr<Map> map_ptr_;
        //std::shared_ptr<PlannerConnection<nx,nu>> planner_connection_;/// connecting segment

        StateType init_state_;///start state
        StateType target_state_;///target state
        double goal_radius_;///goal radius

        double best_cost_ = std::numeric_limits<double>::max();
        std::pair<std::shared_ptr<TreeNodeType>, std::shared_ptr<TreeNodeType>> best_goal_;///best goals
        std::pair<unsigned long, unsigned long> number_of_nodes_;///nodes on kdtree
        bool is_bi_tree_planner_ = false; ///flag to indicate bi-tree
        bool is_optimized_connect_ = false;///flag to indicate adopting BVP to connect two trees
        double optimization_distance_;
        /// a flag to indicate the planner is stopped. This flag is used to terminate long time process
        std::atomic_bool run_flag_;

        std::vector<StateType> connect_states_;
        std::vector<ControlType> connect_control_;
        std::vector<double> connect_duration_;

        void set_opitmization(bool is_optimized,double distance){
            is_optimized_connect_ = is_optimized;
            optimization_distance_ = distance;
        }

        void set_bi_tree(bool bi){
            is_bi_tree_planner_ = bi;
        }

        /***
         * check whether a node is in path
         * @param node
         * @return
         */
        bool is_in_solution_path(std::shared_ptr<TreeNodeType> node) {
            if (best_goal_.first == nullptr) return false;
            if (node->getTreeId() == TreeId::forward) {
                auto it = best_goal_.first;
                while (it) {
                    if (it == node)return true;
                    it = it->getParent();
                }
            } else if (node->getTreeId() == TreeId::backward) {
                auto it = best_goal_.second;
                while (it) {
                    if (it == node)return true;
                    it = it->getParent();
                }
            }
            return false;
        }


public:
/*
        std::vector<StateType> getSolutionVectors() {


            if (best_goal_.first == nullptr) {
                return std::vector<StateType>{};
            }
            auto node = best_goal_.first;
            size_t count = 0;
            while (node) {
                node = node->getParent();
                count++;
            }

            std::vector<StateType> forward_states(count);
            node = best_goal_.first;
            size_t i =1;
            while (node) {
                forward_states[count-i] = node->getState();
                node = node->getParent();
                i++;
            }
            

            std::copy(connect_states_.begin()+1,connect_states_.end()-1,std::back_inserter(forward_states));

            node = best_goal_.second;            
            while (node) {
                forward_states.push_back(node->getState());
                node = node->getParent();
            }
            // std::cout<<count<<"\t"<<forward_states.size()<<std::endl;
            return forward_states;
           
        }

*/

        Planner() = delete;

        /***
         * construct
         * @param system
         */
        explicit Planner(std::shared_ptr<DynamicsModel> dynamics_model_ptr,std::shared_ptr<Map> map_ptr,double goal_radius)
            :dynamics_model_ptr_(dynamics_model_ptr),
            map_ptr_(map_ptr),
            run_flag_(false),
            goal_radius_(goal_radius){
        }

        virtual ~Planner() = default;

        /***
         * necessary setups when starting a planner
         */
        virtual void setup() = 0;
        

        /***
         * one forward step. this function can be called in a while loop for a forward propagating process
         */
        virtual void forward_step() = 0;

        /***
         * one backward step. this function can be called in a while loop for a backward propagating process
         */
        virtual void backward_step() = 0;

        /***
         * one connecting step. this function can be called in a while loop for a connecting process
         */
        virtual void connecting_step() = 0;

        virtual void reset()=0;

        /***
         * get the current numbers of nodes
         * @return a pair of nodes on forward tree and nodes on backward tree
         */
        virtual std::pair<unsigned long, unsigned long> get_node_count() {
            return number_of_nodes_;
        }

        /***
         *
         * @return goal radius
         */
        virtual double get_goal_radius() const {
            return goal_radius_;
        }

        /***
         * set goal radius
         */
        virtual void set_goal_radius(double goal_radius) {
            goal_radius_ = goal_radius;
        }

        /***
         * get start state
         * @return
         */
        virtual StateType get_start_state() const {
            return init_state_;
        }

        /***
         * set start state
         */
        virtual void set_start_state(const StateType &init_state) {
            init_state_ = init_state;
        }

        /***
         * get target state
         * @return
         */
        virtual StateType get_target_state() const {
            return target_state_;
        }

        /***
         * set target state
         */
        virtual void set_target_state(const StateType &target_state) {
            target_state_ = target_state;
        }


        /***
         * get the solution data for solution update observers.
         * @param forward_states
         * @param backward_states
         * @param connect_states
         * @param forward_control
         * @param backward_control
         * @param connect_control
         * @param forward_durations
         * @param backward_durations
         * @param connect_durations
         */
        /*
        std::vector<StateType> getSolutionVectors(double hold_distance = 0.1) {
            std::vector<StateType> forward_states;
            std::vector<StateType> backward_states;
            std::vector<StateType> connect_states;,
            std::vector<ControlType> forward_control;
            std::vector<ControlType> backward_control;
            std::vector<ControlType> connect_control;
            std::vector<double> forward_durations;
            std::vector<double> backward_durations;
            std::vector<double> connect_durations;
           
            if (this->_best_goal.first == nullptr) {
                return;
            }
            auto node = this->_best_goal.first;
            while (node) {
                forward_states.push_back(node->getState());
                forward_control.push_back(node->getEdgeControl());
                forward_durations.push_back(node->getEdgeDuration());
                node = node->getParent();
            }

            node = _best_goal.second;
            backward_durations.push_back(0);
            backward_control.push_back(ControlType());
            while (node) {
                backward_states.push_back(node->getState());
                backward_control.push_back(node->getEdgeControl());
                backward_durations.push_back(node->getEdgeDuration());
                node = node->getParent();
            }
            backward_durations.pop_back();
            backward_control.pop_back();

            if (_planner_connection != nullptr && _planner_connection->_end_states == _best_goal) {
                connect_states = _planner_connection->_states;
                connect_control = _planner_connection->_controls;
                connect_durations = _planner_connection->_durations;
            }

            std::reverse(forward_states.begin(), forward_states.end());
            std::reverse(forward_control.begin(), forward_control.end());
            std::reverse(forward_durations.begin(), forward_durations.end());

            auto distance = [this](const StateType & s1,const StateType& s2){
                return sqrt((s1[0]-s2[0])*(s1[0]-s2[0])+(s1[1]-s2[1])*(s1[1]-s2[1]));
            };


            auto current_state = forward_states.front();
            forward_grid.addVector(current_state*1e6,0.0);

            double t = 0.0;
            for(auto i=0;i<forward_states.size()-1;++i){
                auto temp_state = forward_states[i];
                //std::cout<<"Start State:"<<current_state.transpose()<<'\n';
                auto temp_control=forward_control[i+1];
                int steps = forward_durations[i+1]/_dynamic_system->getStepSize();
                for(auto j=0;j<steps;++j){
                    temp_state = forward(temp_state,temp_control,_dynamic_system->getStepSize());
                    t+=_dynamic_system->getStepSize();
                    if(distance(temp_state,current_state)>hold_distance){
                        current_state = temp_state;
                        forward_grid.addVector(current_state*1e6, forward_grid.getLastTime() + t);
                        t = 0.0;
                    }
                }
            }

            if(!connect_states.empty()) {
                for (auto i = 0; i < connect_states.size() - 1; ++i) {
                    auto temp_state = connect_states[i];
                    auto temp_control = connect_control[i + 1];
                    int steps = connect_durations[i+1] / _dynamic_system->getStepSize();
                    for (auto j = 0; j < steps; ++j) {
                        temp_state = forward(temp_state, temp_control, _dynamic_system->getStepSize());
                        t+=_dynamic_system->getStepSize();
                        if(maxDistance(temp_state,current_state)>SystemConfig::max_distance){
                            current_state = temp_state;
                            forward_grid.addVector(current_state*1e6, forward_grid.getLastTime() + t);
                            t = 0.0;
                        }
                    }
                }
            }


            backward_grid.addVector(backward_states.back(),0.0);
            for(int i=backward_states.size()-1;i>0;--i){
                auto current_state = backward_states[i];
                auto current_control=backward_control[i];
                int steps = backward_durations[i]/_dynamic_system->getStepSize();
                for(auto j=0;j<steps;++j){
                    current_state = forward(current_state,current_control,_dynamic_system->getStepSize());
                    backward_grid.addVector(current_state,backward_grid.getLastTime()+_dynamic_system->getStepSize());
                }
            }
            for(int i=backward_grid.getNumPoints()-1;i>=0;--i){
                t+=_dynamic_system->getStepSize();
                if(maxDistance(backward_grid.getVector(i),current_state)>25e-6){
                    current_state = backward_grid.getVector(i);
                    forward_grid.addVector(current_state*1e6, forward_grid.getLastTime() + t);
                    t = 0.0;
                }
            }

            if(t>_dynamic_system->getStepSize()){
                forward_grid.addVector(backward_grid.getVector(0)*1e6, forward_grid.getLastTime()+t);
            }

            std::stringstream ss;
            forward_grid.print(ss,"","","",10,6,"\t","\n");
            solution_string=ss.str();
        }*/

        /***
         * update best cost
         */
        /*
        virtual void updateBestCost() {
            if (best_goal_.first == nullptr || best_goal_.second == nullptr){
                best_cost_ = std::numeric_limits<double>::max();
            }else {
                double cost = best_goal_.first->getCost() + best_goal_.second->getCost();
                if (planner_connection_ != nullptr && planner_connection_->_end_states.first == best_goal_.first
                    && planner_connection_->_end_states.second == best_goal_.second) {
                    cost += planner_connection_->getTotalCost();
                }
                best_cost_ = cost;
            }
        }*/

        double get_best_cost() const{
            return best_cost_;
        }

        /***
         * notify the planner to stop. here the run_flag can be set to false. Also can call a function to stop the process in dynamic system
         */
        virtual void stop() {
            run_flag_.store(false);
            //dynamic_system_->stop();            
        }

        virtual bool is_running(){
            return run_flag_.load();
        }

        bool is_find_solution(){
            return best_cost_<100000.0;
        }

    };
}


#endif //NANOWIREPLANNER_NANOWIRE_PLANNER_HPP