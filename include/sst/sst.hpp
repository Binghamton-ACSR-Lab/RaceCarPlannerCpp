//
// Created by acsr on 4/27/21.
//

#ifndef ACSR_SST_HPP
#define ACSR_SST_HPP
#include "planner.hpp"
#include <mutex>
#include <unordered_map>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include "array_adaptor.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include "node_adaptor.hpp"


BOOST_GEOMETRY_REGISTER_ARRAY_CS(cs::cartesian)
using nlohmann::json;

namespace acsr{

    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    
    //BOOST_GEOMETRY_REGISTER_NODE_CS(cs::cartesian)

    /***
     * hash function for key = IVector
     */
    template<size_t nc>
    struct GridHash
    {
        std::size_t operator()(const std::array<size_t ,nc> & k) const
        {
            size_t seed = 0;
            for (auto elem:k) {
                seed ^= std::hash<size_t>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    template<typename DynamicsModel,typename Map>
    struct SST : virtual public Planner<DynamicsModel,Map> {
    protected:
    public:
        using Planner<DynamicsModel,Map>::nx_;
        using Planner<DynamicsModel,Map>::nu_;

        constexpr static size_t nc_=DynamicsModel::nc_;
        constexpr static size_t ni_=DynamicsModel::ni_;

        using Planner<DynamicsModel,Map>::dynamics_model_ptr_;
        using Planner<DynamicsModel,Map>::map_ptr_;
        using Planner<DynamicsModel,Map>::best_goal_;
        using Planner<DynamicsModel,Map>::init_state_;
        using Planner<DynamicsModel,Map>::target_state_;
        using Planner<DynamicsModel,Map>::best_cost_;
        using Planner<DynamicsModel,Map>::run_flag_;
        using Planner<DynamicsModel,Map>::optimization_distance_;
        using Planner<DynamicsModel,Map>::connect_states_;
        using Planner<DynamicsModel,Map>::connect_control_;
        using Planner<DynamicsModel,Map>::connect_duration_;

        using StateType = typename Planner<DynamicsModel,Map>::StateType;
        using ControlType = typename Planner<DynamicsModel,Map>::ControlType;

        using NodeType = typename Planner<DynamicsModel,Map>::NodeType;
        using NodePtr = typename Planner<DynamicsModel,Map>::NodePtr;
        using TreeNodeType = typename Planner<DynamicsModel,Map>::TreeNodeType;
        using TreeNodePtr = typename Planner<DynamicsModel,Map>::TreeNodePtr;

        //using Planner<T,STATE_DIMENSION,CONTROL_DIMENSION>::StateType;
        //using Planner<T,STATE_DIMENSION,CONTROL_DIMENSION>::ControlType;
        //using StateType = std::array<T,STATE_DIMENSION>;
        //using ControlType = std::array<T,CONTROL_DIMENSION>;
        using TreeIndexType = std::array<double,ni_>;
        using CellType = std::array<size_t ,nc_>;

        //using NodeType = Node<T,STATE_DIMENSION>;
        //using NodePtr = std::shared_ptr<NodeType>;

        //using TreeNodeType = TreeNode<T,STATE_DIMENSION,CONTROL_DIMENSION>;
        //using TreeNodePtr=std::shared_ptr<TreeNode<T,STATE_DIMENSION,CONTROL_DIMENSION>>;
 
        using SSTNodeType = SSTTreeNode<double,nx_,nu_>;
        using SSTNodePtr = std::shared_ptr<SSTNodeType>;

    protected:
        /***
         * override function
         * @param init_state
         * @param control
         * @param step_length
         * @param steps
         * @param result_state
         * @param duration
         * @return
         */
        virtual std::pair<bool,StateType> forward_propagate(const StateType& init_state, const ControlType& control,double dt, int steps){
            return dynamics_model_ptr_->forward(init_state,control,dt,steps);
        }

        /***
         * override function
         * @param init_state
         * @param control
         * @param step_length
         * @param steps
         * @param result_state
         * @param duration
         * @return
         */
        virtual std::pair<bool,StateType> backward_propagate(const StateType& init_state, const ControlType & control,double dt, int steps){
            return dynamics_model_ptr_->backward(init_state,control,dt,steps);
        }

        /***
         * add a node to tree
         * @param tree_id tree id
         * @param parent the parent node
         * @param state the state of the node
         * @param control the contrl of the edge
         * @param duration the duration of the edge
         */
        virtual SSTNodePtr add_to_tree(TreeId tree_id,
                                     TreeNodePtr parent,
                                     const StateType& state,
                                     const ControlType& control,
                                     double duration){
            auto target_node = tree_id==TreeId::forward?goal_:root_;
            if(best_cost_ <= parent->get_cost() + duration + dynamics_model_ptr_->getHeuristic(state,target_node->get_state()))
                return nullptr;

            
            //std::cout<<"Cell Dimension "<<CELL_DIMENSION<<std::endl;
            auto& map = (tree_id==TreeId::forward?forward_prox_map:backward_prox_map);
            CellType iv = dynamics_model_ptr_->state_to_cell(map_ptr_,state);
            
            
            SSTNodePtr previous_node = nullptr;
            if(map.find(iv)!=map.end() ) {
                if(!(map[iv].expired()))
                    previous_node = map[iv].lock();
            }
            
            if(previous_node
               && previous_node->get_cost() <= parent->get_cost() + duration)
                return nullptr;
            
            ///create a new tree node
            auto new_node = std::make_shared<SSTNodeType>(tree_id,state);
            new_node->set_edge({control,duration});
            new_node->set_parent(parent);
            
            new_node->set_cost(parent->get_cost() + duration);
            parent->add_child(new_node);
            
           
            /// the original witness has a monitor node but that node has greater cost than the node we are working
            if(previous_node){
                previous_node->set_active(false);
                removePointFromSearchTree(previous_node);
                auto best_goal = tree_id==TreeId::forward?best_goal_.first:best_goal_.second;
                if(previous_node!=best_goal) {
                    while (previous_node->get_children().empty() && !previous_node->is_active()) {
                        auto next = previous_node->get_parent();
                        removeLeaf(previous_node);
                        previous_node = std::static_pointer_cast<SSTNodeType>(
                                next);
                    }
                }
            }



            ///set new node as the current active point
            map[iv] = new_node;
            new_node->set_active(true);
            add_node_to_search_tree(new_node);
            
            return new_node;
        }

        /***
         * check the distance of a new node to another tree
         * @param node
         */
        void check_connection(const TreeNodePtr& node){
            if(node == nullptr)return;
            
            TreeId another_treeid = (node->getTreeId()==TreeId::forward?TreeId::backward:TreeId::forward);

            auto nearest_vec_node = get_near_node_by_radius(node->get_state(),another_treeid,optimization_distance_);
            if(!nearest_vec_node.empty()){
                for(auto& n : nearest_vec_node){
                    auto pair_node = std::static_pointer_cast<TreeNodeType>(n);
                    if(pair_node->get_cost()+node->get_cost() >best_cost_){
                        continue;
                    }
                    ///distance within goal radius, trying to update solution
                    auto d = this->dynamics_model_ptr_->euclidean_distance(pair_node->get_state(),node->get_state());
                    if(d<this->goal_radius_){
                        if(node->getTreeId() == TreeId::forward) {
                            this->best_goal_.first = node;
                            this->best_goal_.second = std::static_pointer_cast<TreeNodeType>(n);
                        }else{
                            this->best_goal_.first = std::static_pointer_cast<TreeNodeType>(n);
                            this->best_goal_.second = node;
                        }
                        best_cost_ = best_goal_.first->get_cost() + best_goal_.second->get_cost();
                        //branchBound(root_);
                        //branchBound(goal_);
                    }else{
                        if(node->getTreeId()==TreeId::forward) {
                            this->optimize_set[0] = std::make_pair(node,pair_node);
                        }else{
                            this->optimize_set[1] = std::make_pair(pair_node,node);
                        }
                    }
                }

            }
        }

        /***
         * add node to kdtree
         * @param node
         */
        virtual void add_node_to_search_tree(TreeNodePtr node){
            if(node->getTreeId()==TreeId::forward){
                std::scoped_lock<std::mutex> lock1(forward_tree_mutex);
                //backward_rtree.insert({node->get_state(),node});
                forward_rtree.insert({dynamics_model_ptr_->state_to_index(node->get_state()),node});
                this->number_of_nodes_.first+=1;
            }
            else if(node->getTreeId()==TreeId::backward){
                std::scoped_lock<std::mutex> lock1(backward_tree_mutex);
                //backward_rtree.insert({node->get_state(),node});
                backward_rtree.insert({dynamics_model_ptr_->state_to_index(node->get_state()),node});
                this->number_of_nodes_.second+=1;
            }

        }

        /***
         * remove node from tree
         * @param node
         */
        virtual void remove_node_from_search_tree(TreeNodePtr node){
            if(node== nullptr )
                return;

            if(node->getTreeId()==TreeId::forward){
                std::scoped_lock<std::mutex> lock1(forward_tree_mutex);
                auto p = forward_rtree.remove({dynamics_model_ptr_->state_to_index(node->get_state()),node});
                this->number_of_nodes_.first -= p;
            }
            else if(node->getTreeId()==TreeId::backward){
                std::scoped_lock<std::mutex> lock1(backward_tree_mutex);
                auto p = backward_rtree.remove({dynamics_model_ptr_->state_to_index(node->get_state()),node});
                this->number_of_nodes_.second -= p;
            }
        }

        /***
         * get nearest n nodes
         * @param state
         * @param tree_id
         * @param count
         * @return
         */
        virtual std::vector<NodePtr> get_near_node_by_count(const StateType & state,TreeId tree_id,int count){
            std::vector<std::pair<TreeIndexType,NodePtr>> temp_vec;
            
            if(tree_id == TreeId::forward){
                std::scoped_lock<std::mutex> lock(forward_tree_mutex);
                forward_rtree.template query(bgi::nearest(dynamics_model_ptr_->state_to_index(state),count),std::back_inserter(temp_vec));
                //forward_rtree.template query(bgi::nearest(state,count),std::back_inserter(temp_vec));
            }else if(tree_id == TreeId::backward){
                std::scoped_lock<std::mutex> lock(backward_tree_mutex);
                backward_rtree.template query(bgi::nearest(dynamics_model_ptr_->state_to_index(state),count),std::back_inserter(temp_vec));
            }
            std::vector<NodePtr> vec(temp_vec.size());
            std::transform(temp_vec.begin(),temp_vec.end(),vec.begin(),[](auto& t){
                return t.second;
            });

            return vec;
        }

        /***
         * get node within a radius arount state
         * @param state
         * @param tree_id
         * @param radius
         * @return first: the vector containing all nodes within radius around state, second: the nearest node, no matter whether it's in radius
         */
        virtual std::vector<NodePtr> get_near_node_by_radius(const StateType& state,TreeId tree_id,double radius){
            std::vector<NodePtr> vec;
            if(tree_id == TreeId::forward){
                std::scoped_lock<std::mutex> lock(forward_tree_mutex);
                /*forward_rtree.template query(bgi::nearest(state,20) && bgi::satisfies([&](const NodePtr& node){
                    return bg::distance(node,state)<radius;
                }),std::back_inserter(vec));*/
               // for ( auto it = forward_rtree.qbegin(bgi::nearest(state, 20)) ;
                for ( auto it = forward_rtree.qbegin(bgi::nearest(dynamics_model_ptr_->state_to_index(state), 20)) ;
                      it != forward_rtree.qend() ; ++it )
                {
                    auto n = std::static_pointer_cast<SSTNodeType>(it->second);
                    if(n == nullptr){
                        forward_rtree.remove(*it);
                    }
                    else if(!(n->isActive()))forward_rtree.remove(*it);
                    else{
                        if(dynamics_model_ptr_->euclideanDistance(n->get_state(),state)<radius)
                            vec.push_back(it->second);
                            // std::cout<<"it value is:"<<vec<<std::endl;
                    }
                }
            }else if(tree_id == TreeId::backward){
                std::scoped_lock<std::mutex> lock(backward_tree_mutex);
                /*
                backward_rtree.template query(bgi::nearest(state,20) && bgi::satisfies([&](const NodePtr& node){
                    if(node == nullptr)return false;
                    if(node->get_state().data() == nullptr)return false;
                    return (node->get_state()-state).norm()<radius;
                }),std::back_inserter(vec));*/
                for ( auto it = backward_rtree.qbegin(bgi::nearest(dynamics_model_ptr_->state_to_index(state), 20)) ;
                      it != backward_rtree.qend() ; ++it )
                {
                    auto n = std::static_pointer_cast<SSTNodeType>(it->second);
                    if(n== nullptr){
                        backward_rtree.remove(*it);
                    }
                    else if(!(n->isActive()))backward_rtree.remove(*it);
                    else{
                        if(dynamics_model_ptr_->euclideanDistance(n->get_state(),state)<radius)
                            vec.push_back(it->second);
                    }
                }
            }
            return vec;
        }

        /***
         * select a pair node for optimizing connecting
         * @param node
         * @param pair_node
         * @return
         */
        std::pair<double,TreeNodePtr> choose_pair_node(TreeNodePtr node){
            double min_value=-1.0;
            if(node== nullptr )return -1;
            TreeId other_tree_id = node->getTreeId()==TreeId::backward?TreeId::forward:TreeId::backward;
            auto nearest_vec = getNearNodeByRadius(node->get_state(),other_tree_id,optimization_distance_);
            if(nearest_vec.empty()){
                return {min_value, nullptr};
            }

            TreeNodePtr pair_node{nullptr};
            min_value = std::numeric_limits<double>::max();
            for(auto temp_pair_node:nearest_vec){
                auto temp_tree_node = std::dynamic_pointer_cast<SSTNodeType>(temp_pair_node);
                if(!temp_tree_node->isActive())
                    continue;
                auto temp_cost = this->dynamics_model_ptr_->getHeuristic(temp_tree_node->get_state(),node->get_state()) + temp_tree_node->get_cost();
                if( temp_cost < min_value){
                    min_value=temp_cost;
                    pair_node = temp_tree_node;
                }
            }
            if(pair_node== nullptr)return {-1,nullptr};
            else return {min_value + node->get_cost()+dynamics_model_ptr_->get_heuristic(node->get_state(),pair_node->get_state()),pair_node};
        }

        /***
         * remove leaf node
         * @param node
         */
        virtual void remove_leaf(TreeNodePtr node){
            if(node->get_parent()){
                node->get_parent()->remove_child(node);
            }
        }

        /***
         * branch trim
         * @param node
         */
        virtual void branch_bound(TreeNodePtr node){
            std::lock_guard<std::recursive_mutex> lk(branch_mutex);
            auto children = node->get_children();
            auto temp_target = node->get_tree_id()==TreeId::forward?target_state_:init_state_;
            if(node->get_cost() + dynamics_model_ptr_->get_heuristic(node->get_state(),temp_target)  > best_cost_){
                remove_branch(node);
            }else{
                for(auto &n:children){
                    branch_bound(n);
                }
            }
        }

        /***
         * remove a branch
         * @param node the branch root going to be removed
         */
        virtual void remove_branch(TreeNodePtr node){
            auto children = node->get_children();
            if(node->get_parent()!=nullptr){
                node->get_parent()->remove_child(node);
            }
            auto sst_node = std::static_pointer_cast<SSTNodeType>(node);
            if(sst_node->is_active()) {
                sst_node->set_active(false);
                remove_node_from_search_tree(node);
            }
            for(auto& n:children){
                remove_branch(n);
            }
        }

    protected:
        bgi::rtree< std::pair<TreeIndexType,NodePtr>, bgi::quadratic<32> > forward_rtree;
        bgi::rtree< std::pair<TreeIndexType,NodePtr>, bgi::quadratic<32> > backward_rtree;
        //bgi::rtree< NodePtr, bgi::quadratic<32> > forward_rtree;
        //bgi::rtree< NodePtr, bgi::quadratic<32> > backward_rtree;

        std::mutex forward_tree_mutex;
        std::mutex backward_tree_mutex;

        std::mutex optimize_set_mutex;
        std::map<double,std::pair<std::weak_ptr<TreeNodeType>,std::weak_ptr<TreeNodeType>>> optimize_set;

        SSTNodePtr root_;
        SSTNodePtr goal_;

        std::unordered_map<CellType,std::weak_ptr<SSTNodeType>,GridHash<nc_>> forward_prox_map; ///map storing forward active node but in weak_ptr
        std::unordered_map<CellType,std::weak_ptr<SSTNodeType>,GridHash<nc_>> backward_prox_map; ///map storing backword active node but in weak_ptr

        std::recursive_mutex branch_mutex;

        double step_size_,sst_delta_near_;
        int min_time_steps_,max_time_steps_;
        
    public:
        SST() = delete;

        /***
         * constructor
         * @param dynamic_system
         */
        explicit SST(std::shared_ptr<DynamicsModel> dynamic_model_ptr,std::shared_ptr<Map> map_ptr,double goal_radius, const json& sst_params):
            Planner<DynamicsModel,Map>(dynamic_model_ptr,map_ptr,goal_radius)
        {
            step_size_ = sst_params.at("step_size");
            min_time_steps_ = sst_params.at("min_time_steps");
            max_time_steps_ = sst_params.at("max_time_steps");
            sst_delta_near_ = sst_params.at("sst_delta_near");
            if(sst_params.contains("optimization_distance"))
                optimization_distance_ = sst_params.at("optimization_distance");
            else
                optimization_distance_=goal_radius;
        }

        /***
         * destructor
         */
        virtual ~SST() override{
            forward_rtree.clear();
            forward_rtree.clear();
            optimize_set.clear();
            forward_prox_map.clear();
            backward_prox_map.clear();
            root_ = nullptr;
            goal_ = nullptr;
        }

        virtual void reset() override{
            forward_rtree.clear();
            forward_rtree.clear();
            optimize_set.clear();
            forward_prox_map.clear();
            backward_prox_map.clear();
            root_ = nullptr;
            goal_ = nullptr;
            this->number_of_nodes_=std::make_pair(0,0);
            this->best_goal_ = std::make_pair(nullptr,nullptr);
            this->best_cost_ = std::numeric_limits<double>::max();
            this->run_flag_ = false;
        }

        /***
         * override setup
         */
        virtual void setup() override{
            run_flag_ = true;

            ///init root and goal
            root_ = std::make_shared<SSTNodeType>(TreeId::forward, init_state_);
            root_->set_edge(TreeEdge<double,nu_>(ControlType{},0.0));
            root_->set_active(true);
            add_node_to_search_tree(root_);
            forward_prox_map[dynamics_model_ptr_->state_to_cell(map_ptr_, init_state_)]= root_;

            goal_ = std::make_shared<SSTNodeType>(TreeId::backward, target_state_);
            goal_->set_edge(TreeEdge<double,nu_>(ControlType{},0.0));
            goal_->set_active(true);
            add_node_to_search_tree(goal_);
            backward_prox_map[dynamics_model_ptr_->state_to_cell(map_ptr_,target_state_)]= goal_;
        }

        /***
         * override forward step
         */
        virtual void forward_step() override{
            auto temp_state = dynamics_model_ptr_->random_state(map_ptr_);
            auto temp_control = dynamics_model_ptr_->random_control();
            auto nearest_vec_node = get_near_node_by_radius(temp_state,TreeId::forward,sst_delta_near_);
            NodePtr parent;
            
            if(!nearest_vec_node.empty()) {
                parent = *std::min_element(nearest_vec_node.begin(), nearest_vec_node.end(),
                                           [](const NodePtr &node1, const NodePtr &node2) 
                                           {
                                               return std::static_pointer_cast<TreeNodeType>(node1)->get_cost() < std::static_pointer_cast<TreeNodeType>(node2)->get_cost();
                                           });
            }else{
                parent = get_near_node_by_count(temp_state,TreeId::forward,1).front();
            }

            auto steps = acsr::Random::get_instance().random_value(min_time_steps_,max_time_steps_);
            auto result = forward_propagate(parent->get_state(),temp_control,step_size_,steps);
            //double duration;            
            if(result.first){
                auto duration = step_size_*steps;
                auto new_node = add_to_tree(TreeId::forward,std::static_pointer_cast<TreeNodeType>(parent),result,temp_control,duration);
                check_connection(new_node);
            }
        }

        /***
         * override backward step
         */
        virtual void backward_step() override{

            auto temp_state = dynamics_model_ptr_->random_state(map_ptr_);
            auto temp_control = dynamics_model_ptr_->random_control();
            auto nearest_vec_node = get_near_node_by_radius(temp_state,TreeId::backward,sst_delta_near_);
            NodePtr parent;
            
            if(!nearest_vec_node.empty()) {
                parent = *std::min_element(nearest_vec_node.begin(), nearest_vec_node.end(),
                                           [](const NodePtr &node1,const NodePtr &node2) {
                                               return std::static_pointer_cast<TreeNodeType>(node1)->get_cost() < std::static_pointer_cast<TreeNodeType>(node2)->get_cost();
                                           });
            }else{
                parent = get_near_node_by_count(temp_state,TreeId::backward,1).front();
            }
            auto steps = acsr::Random::random_value(min_time_steps_,max_time_steps_);
            auto result=backward_propagate(parent->get_state(),temp_control,step_size_,steps);
            //double duration;
            if(result.first){
                auto duration = step_size_*steps;
                auto new_node = add_to_tree(TreeId::backward,std::static_pointer_cast<TreeNodeType>(parent),result,temp_control,duration);
                check_connection(new_node);
            }
        }

        /***
         * override connecting step
         */
        void connecting_step() override{
            /*
            if(!run_flag_)return;            
            if(this->optimize_set.empty())
                return;
            TreeNodePtr explore_node = nullptr;
            TreeNodePtr temp_explore_node = nullptr;
            TreeNodePtr target = nullptr;
            TreeNodePtr temp_target = nullptr;

            auto it = this->optimize_set.begin();
            while(it!=this->optimize_set.end()){
                if(!this->run_flag_)return;
                std::scoped_lock<std::mutex> lock(optimize_set_mutex);
                if(it->second.first.expired() || it->second.second.expired()){
                    it = optimize_set.erase(it);
                    continue;
                }
                explore_node = it->second.first.lock();
                target = it->second.second.lock();
                optimize_set.erase(it);
                break;
            }
            auto h1 = chooseOtherNearestForOptimization(explore_node,temp_target);
            auto h2 = chooseOtherNearestForOptimization(temp_target,temp_explore_node);
            if(h1<0 && h2<0)return;
            double total_cost;
            if(h1<0 && h2>0){
                explore_node=temp_explore_node;
                total_cost=h2;
            }else if(h1>0 && h2<0){
                target=temp_target;
                total_cost=h1;
            }else if(h1>0 && h1<h2){
                target=temp_target;
                total_cost=h1;
            }else{
                explore_node=temp_explore_node;
                total_cost=h2;
            }

            if(explore_node ==nullptr
                || !std::static_pointer_cast<SSTNodeType>(explore_node)->isActive()
                || target==nullptr
                || !std::static_pointer_cast<SSTNodeType>(target)->isActive()
                || total_cost > best_cost_){
                return;
            }

            std::vector<StateType> vec_state;
            std::vector<ControlType> vec_control;
            std::vector<double> vec_duration;

            auto init_state = explore_node->get_state();
            auto target_state = target->get_state();

            bool optimized = dynamics_model_ptr_->connect(init_state, target_state,
                                                            vec_state, vec_control,vec_duration);
            if(!run_flag_)return;
            if(explore_node == nullptr || target == nullptr)
                return;

            if(!optimized){
                return;
            }
            if(vec_duration.size()==0)
                return;

            auto total_duration = std::accumulate(vec_duration.begin(),vec_duration.end(),0.0);


            if(explore_node->get_cost() + target->get_cost() + total_duration < best_cost_){

                std::cout<<"Explore new solution\n";
                best_goal_.first = explore_node;
                best_goal_.second = target;
                best_cost_ = explore_node->get_cost() + target->get_cost() + total_duration;
                connect_states_ = vec_state;
                connect_duration_ = vec_duration;
                connect_control_ = vec_control;
                //this->updateBestCost();
                branchBound(root_);
                branchBound(goal_);
                
            }*/
        }

        std::vector<StateType> get_solution_vectors() {

            std::vector<StateType> forward_states;
            std::vector<StateType> backward_states;
            std::vector<ControlType> forward_control;
            std::vector<ControlType> backward_control;
            std::vector<double> forward_durations;
            std::vector<double> backward_durations;

            if (this->best_goal_.first == nullptr) {
                return std::vector<StateType>{};
            }
            auto node = this->best_goal_.first;
            while (node) {
                forward_states.push_back(node->get_state());
                forward_control.push_back(node->get_edge_control());
                forward_durations.push_back(node->get_edge_duration());
                node = node->get_parent();
            }


            std::reverse(forward_states.begin(), forward_states.end());
            std::reverse(forward_control.begin(), forward_control.end());
            std::reverse(forward_durations.begin(), forward_durations.end());

            std::vector<StateType> forward_result_states{forward_states.front()};
            for(auto i=0;i<forward_states.size()-1;++i){
                auto temp_state = forward_states[i];
                auto temp_control=forward_control[i+1];
                int steps = forward_durations[i+1]/step_size_;
                for(auto j=0;j<steps;++j){
                    auto result = dynamics_model_ptr_->forward(temp_state,temp_control,step_size_,1);
                    forward_result_states.push_back(result.second);
                }
            }

            node = this->best_goal_.second;
            backward_durations.push_back(0.0);
            backward_control.push_back(ControlType());
            while (node) {
                backward_states.push_back(node->get_state());
                backward_control.push_back(node->get_edge_control());
                backward_durations.push_back(node->get_edge_duration());
                node = node->getParent();
            }

            std::reverse(backward_states.begin(), backward_states.end());
            std::reverse(backward_control.begin(), backward_control.end());
            std::reverse(backward_durations.begin(), backward_durations.end());
            std::vector<StateType> backward_result_states{backward_states.front()};

            for(auto i=0;i<backward_states.size()-1;++i){
                auto temp_state = backward_states[i];
                auto temp_control=backward_control[i+1];
                int steps = backward_durations[i+1]/step_size_;
                for(auto j=0;j<steps;++j){
                    auto result = dynamics_model_ptr_->backward(temp_state,temp_control,step_size_,1);
                    backward_result_states.push_back(result.second);
                }
            }
            forward_result_states.insert( forward_result_states.end(),
                                          std::make_move_iterator(backward_result_states.rbegin()),
                                          std::make_move_iterator(backward_result_states.rend() ));
            return forward_result_states;

        }


        std::vector<StateType> get_solution_vectors(double hold_distance) {

            std::vector<StateType> forward_states;
            std::vector<StateType> backward_states;
            std::vector<ControlType> forward_control;
            std::vector<ControlType> backward_control;
            std::vector<double> forward_durations;
            std::vector<double> backward_durations;
           
            if (this->best_goal_.first == nullptr) {
                return std::vector<StateType>{};
            }
            auto node = this->best_goal_.first;
            while (node) {
                forward_states.push_back(node->get_state());
                forward_control.push_back(node->get_edge_control());
                forward_durations.push_back(node->get_edge_duration());
                node = node->get_parent();
            }


            std::reverse(forward_states.begin(), forward_states.end());
            std::reverse(forward_control.begin(), forward_control.end());
            std::reverse(forward_durations.begin(), forward_durations.end());                       
            std::vector<StateType> forward_result_states{forward_states.front()};
            //double t = 0.0;
            double s =0.0;
            for(auto i=0;i<forward_states.size()-1;++i){
                auto temp_state = forward_states[i];
                auto temp_control=forward_control[i+1];
                int steps = forward_durations[i+1]/step_size_;
                for(auto j=0;j<steps;++j){
                    auto result = dynamics_model_ptr_->forward(temp_state,temp_control,step_size_,1);
                    auto& result_state = result.second;
                    s+=dynamics_model_ptr_->euclidean_distance(temp_state,result_state);
                    temp_state = result_state;
                    if(s>hold_distance){
                        forward_result_states.push_back(result_state);
                        s=0.0;
                    }
                }
            }
            if(dynamics_model_ptr_->euclidean_distance(forward_result_states.back(),forward_states.back())<hold_distance/2){
                forward_result_states.push_back(forward_states.back());
            }

            node = this->best_goal_.second;
            backward_durations.push_back(0.0);
            backward_control.push_back(ControlType());
            while (node) {
                backward_states.push_back(node->get_state());
                backward_control.push_back(node->get_edge_control());
                backward_durations.push_back(node->get_edge_duration());
                node = node->getParent();
            }

            std::reverse(backward_states.begin(), backward_states.end());
            std::reverse(backward_control.begin(), backward_control.end());
            std::reverse(backward_durations.begin(), backward_durations.end());
            std::vector<StateType> backward_result_states{backward_states.front()};
            s =0.0;
            for(auto i=0;i<backward_states.size()-1;++i){
                auto temp_state = backward_states[i];
                auto temp_control=backward_control[i+1];
                int steps = backward_durations[i+1]/step_size_;
                std::cout<<"steps: "<<steps<<std::endl;
                for(auto j=0;j<steps;++j){
                    auto result = dynamics_model_ptr_->backward(temp_state,temp_control,step_size_,1);
                    auto& result_state = result.second;
                    s+=dynamics_model_ptr_->euclidean_distance(temp_state,result_state);
                    temp_state = result_state;
                    if(s>hold_distance){
                        backward_result_states.push_back(result_state);
                        s=0.0;
                    }
                }
            }
            if(dynamics_model_ptr_->euclidean_distance(backward_result_states.back(),backward_states.back())<hold_distance/2){
                backward_result_states.push_back(backward_states.back());
            }
            forward_result_states.insert( forward_result_states.end(), backward_result_states.rbegin(), backward_result_states.rend() );
            return forward_result_states;
                  
        }
    


    };
}


#endif //NANOWIREPLANNER_NANOWIRE_SST_HPP
