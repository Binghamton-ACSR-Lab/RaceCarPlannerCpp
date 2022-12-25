//
// Created by acsr on 4/26/21.
//

#ifndef ACSR_PLANNER_TREE_NODE_HPP
#define ACSR_PLANNER_TREE_NODE_HPP

#include <list>
#include <unordered_map>
#include <array>
#include <memory>
#include <stdexcept>
//#include <eigen3/Eigen/Dense>

namespace acsr {
/***
 * denote in which tree a node is
 */
    enum TreeId {
        forward = 0,
        backward,
        connection
    };

    template <class T,int nu>
    using TreeEdge =std::pair<std::array<T,nu>,T>;

    template <class T,int nx,int nu>
    class SSTTreeNode;

    template <class T,int nx,int nu>
    class TreeNode;

    template <class T,int nx>
    class Node;


    template <class T,int nx>
    class Node{
    protected:
        std::array<T,nx> state_;
    public:
        /***
         * delete defaut constructor
         */
        Node() = delete;

        
        explicit Node(const std::array<T,nx>& state):state_(std::move(state)){
        }

        /***
         * copy constructor
         * @param _node
         */
        Node(const Node& node){
            state_ = node.state_;
        }

        /***
         * operator =
         * @param _node
         * @return
         */
        Node& operator=(const Node& node){
            if(&node == this)return *this;
            state_ = node.state_;
            return *this;
        }

        /***
         * override operator []
         * @param index
         * @return
         */
        const T operator[](size_t index) const{
            if(index>nx)
                throw(std::out_of_range("index exceeds state range."));
            return state_[index];
        }

        /***
         * operator []
         * @param index
         * @return
         */
        T& operator[](size_t index) {
            if(index>nx)
                throw(std::out_of_range("index exceeds state range."));
            return state_[index];
        }

        /***
         * destructor
         */
        virtual ~Node()=default;

        /***
         * get the state
         * @return
         */
        virtual std::array<T,nx> get_state() const{
            return state_;
        }

    };

    template <class T,int nx,int nu>
    class TreeNode: public Node<T,nx> {
        //using TreeNodePtr = std::shared_ptr <TreeNode<T,nx,nu>>;
    protected:
        std::weak_ptr <TreeNode> parent_;
        std::list<std::shared_ptr <TreeNode<T,nx,nu>>> children_;
        TreeEdge<T,nu>  edge_;
        T cost_;
        TreeId tree_id_;
        //TreeNodeState _tree_node_state;


    public:
        /***
         * delete default constructor
         */
        TreeNode() = delete;

        /***
         * delete copy constructor
         */
        TreeNode(const TreeNode &) = delete;

        /***
         * delete assign operator
         * @return
         */
        TreeNode &operator=(const TreeNode &) = delete;

        /***
         * constructor
         * @param id: tree id
         * @param pt: node state
         */
        explicit TreeNode(TreeId id, const std::array<T,nx> &pt):Node<T,nx>(pt),tree_id_(id),cost_(0.0) {

        }

        /***
         * deconstructor
         */
        ~TreeNode() override {
            children_.clear();
        }

        /***
         * set the parent of this node
         * @param parent
         */
        void set_parent(const std::shared_ptr <TreeNode<T,nx,nu>> &parent) {
            parent_ = parent;
        }

        /***
         * get parent
         * @return
         */
        std::shared_ptr <TreeNode<T,nx,nu>> get_parent() const {
            if (!parent_.expired())
                return parent_.lock();
            else
                return nullptr;
        }

        /***
         * set node  cost
         */
        void set_cost(T cost) {
            cost_ = cost;
        }

        /***
         * get node cost
         * @return
         */
        T get_cost() const {
            return cost_;
        }

        /***
         * add a child
         * @param child a child node
         */
        void add_child(const std::shared_ptr <TreeNode<T,nx,nu>> &child) {
            children_.emplace_back(child);
        }

        /***
         * remove a child
         * @param child
         */
        void remove_child(std::shared_ptr <TreeNode<T,nx,nu>> &child) {
            children_.remove(child);
        }

        /***
         * clean the children
         */
        void clear_children() {
            children_.clear();
        }

        /***
         * get the children vector
         * @return
         */
        std::list <std::shared_ptr <TreeNode<T,nx,nu>>> get_children() {
            return children_;
        }

        /***
         * get the tree id
         * @return
         */
        TreeId get_tree_id() const {
            return tree_id_;
        }

        /***
         * set tree edge
         * @param _edge
         */
        void set_edge(const TreeEdge<T,nu> &edge) {
            edge_ = edge;
        }

        /***
         * get the control of the edge
         * @return
         */
        std::array<T,nu> get_edge_control() const {
            return edge_.first;
        }

        /***
         * get the duration of the edge
         * @return
         */
        T get_edge_duration() const {
            return edge_.second;
        }

        bool operator==(const TreeNode<T,nx,nu> &node) const {
            return this->state_ == node.state_ && edge_ == node.edge_ && tree_id_ == node.tree_id_ &&
                   cost_ == node.cost_;
        }

    };

    template <class T,int nx,int nu>
    class SSTTreeNode : public TreeNode<T,nx,nu>{

    protected:
        bool active_state_;

    public:

        SSTTreeNode(/* args */) = delete;
        SSTTreeNode(const SSTTreeNode&) = delete;
        SSTTreeNode& operator=(const SSTTreeNode&) = delete;

        /***
         * inherits from treenode
         * @param id
         * @param pt
         */
        SSTTreeNode(TreeId id,const std::array<T,nx>& pt): TreeNode<T,nx,nu>(id,pt){

        }

        /***
         * check whether this sstnode active
         * @return
         */
        bool is_active(){
            return active_state_;
        }

        /***
         * set activeness
         * @param state
         */
        void set_active(bool state){
            active_state_ = state;
        }

        /***
         * desctructor
         */
        ~SSTTreeNode() override =default;
    };


}

#endif //NANOWIREPLANNER_TREE_NODE_HPP