//
// Created by lab on 9/2/22.
//
#include <gtest/gtest.h>
#include "csv_reader.hpp"
#include <casadi/casadi.hpp>
#include "rtree.hpp"

TEST(CSVReaderTest,Reader){
    acsr::CSVReader reader("tracks/temp.csv");
    auto data = reader.read();
    casadi::DM dm;


    EXPECT_TRUE(data.toDM(dm));
    EXPECT_TRUE(dm.columns()==2);
    EXPECT_NEAR(2626.03,double(dm(0,0)),1e-3);
}

TEST(RTreaTest,RTree){
    std::vector<std::vector<double>> pts{{0.0,0.0},{1.0,0.0},{0.0,1.0},{1.0,1.0},{0.5,1.5}};
    acsr::RTree rTree(pts);

    auto [pt,index] = rTree.findNearest({0.8,0.2});

    EXPECT_TRUE(index==1);
    EXPECT_NEAR(pt[0],1.0,1e-3);
    EXPECT_NEAR(pt[1],0.0,1e-3);
}