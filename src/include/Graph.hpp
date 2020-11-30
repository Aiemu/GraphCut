#ifndef GRAPHCUT_GRAPH_HPP
#define GRAPHCUT_GRAPH_HPP

#include <cstdio>
#include <iostream>
#include <queue>
#include <string>

#include <opencv2/opencv.hpp>

// define graph class
class Graph{
public:
    Graph(int, int, int);
    Graph(std::string);
    ~Graph();

    cv::Mat img;
    int get_h();
    int get_w();

private:
    int w;
    int h;
};

Graph::Graph(int h, int w, int i_type){
    // read in img
    img = cv::Mat::zeros(h, w, i_type);

    this->w = w;
    this->h = h;
}

Graph::Graph(std::string path) {
    // read in img
    img = cv::imread(path);

    if (img.empty()) {
        cv::Mat image;
        cv::VideoCapture video_capture;

        image = video_capture.open(path);

        if (!video_capture.isOpened()) {
            video_capture.release();
            return;
        }

        for (; video_capture.read(image); ) {
            img = image.clone();
        }

        if (img.empty()) {
            return;
        }

        video_capture.release();
    }
    w = img.cols;
    h = img.rows;
}

Graph::~Graph() {
    // destroy img
    img.release();
}

int Graph::get_h() {
    return h;
}

int Graph::get_w() {
    return w;
}

#endif //GRAPHCUT_GRAPH_HPP
