//
// Created by acsr on 9/22/22.
//

#ifndef RACECARPLANNER_SHAPE_HPP
#define RACECARPLANNER_SHAPE_HPP

namespace acsr{
    template<typename T>
    class Shape{
    public:
        virtual bool collision(T x,T y)=0;
        virtual T area() = 0;
        virtual std::string type()=0;
        virtual std::pair<std::vector<T>,std::vector<T>> get_plot_data()=0;
        virtual ~Shape()=default;

    };

    template<typename T>
    class Rectangle : public Shape<T>{
    public:
        Rectangle()=delete;
        Rectangle(T x,T y,T width,T height):x_{x},y_{y},width_{width},height_{height}{
        }

        bool collision(T x,T y) override{
            return (x>=x_ && y>=y_ && x<=x_+width_ && y<=y_+height_);
        }

        T area() override{
            return width_*height_;
        }

        std::string type() override{
            return "rectangle";
        }

        std::pair<std::vector<T>,std::vector<T>> get_plot_data() override{
            return std::make_pair<std::vector<T>,std::vector<T>>(std::vector<T>{x_,x_+width_,x_+width_,x_,x_},std::vector<T>{y_,y_,y_+height_,y_+height_,y_});
        }

    private:
        T x_{};
        T y_{};
        T width_{};
        T height_{};
    };

    template<typename T>
    class Circle : public Shape<T>{
    public:
        Circle()=delete;
        Circle(T x,T y,T radius):x_{x},y_{y},radius_{radius}{
        }

        bool collision(T x,T y) override{
            return (x-x_)*(x-x_)+(y-y_)*(y-y_)<=radius_*radius_;
        }

        T area() override{
            return M_PI*radius_*radius_;
        }

        std::string type() override{
            return "circle";
        }

        std::pair<std::vector<T>,std::vector<T>> get_plot_data() override{
            std::vector<T> phi(101);
            std::generate(phi.begin()+1,phi.end(),[n=0.0]() mutable {
                return n+=2*M_PI/100;
            });
            std::vector<T> x(101),y(101);
            std::transform(phi.begin(),phi.end(),x.begin(),[this](T v){
                return x_+radius_*cos(v);
            });
            std::transform(phi.begin(),phi.end(),y.begin(),[this](T v){
                return y_+radius_*sin(v);
            });

            return std::make_pair(x,y);
        }

    private:
        T x_{};
        T y_{};
        T radius_{};
    };





}


#endif //RACECARPLANNER_SHAPE_HPP
