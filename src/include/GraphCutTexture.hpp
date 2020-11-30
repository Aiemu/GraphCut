#ifndef GRAPHCUT_GRAPHCUTTEXTURE_HPP
#define GRAPHCUT_GRAPHCUTTEXTURE_HPP

#define RANDOM 0
#define GLOBAL_BEST 1
#define LOCAL_BEST 2

#define ITERATION_TIME 70
#define MIN_DELTA 0.1
#define HALF_MIN_DELTA 0.05
#define SCALE_SIZE 1000

#include <unistd.h>
#include <iostream>
#include <random>
#include <string>
#include <cstring>
#include <map>
#include <queue>
#include <cmath>

#include <Graph.hpp>
#include <MaxFlow.hpp>
#include <MinCut.hpp>

// GraphCut
class GraphCutTexture {
public:
    GraphCutTexture(std::string path);
    ~GraphCutTexture();

    // acquire private param
    int getModel();
    int getWidth();
    int getHeight();
    void setModel(int model);

    // main process
    bool graphCut(int h, int w, int i_time);

    // output
    void saveTexture(const std::string path);
    void saveSeams(const std::string path_s);

private:
    // compute model
    int model = RANDOM;

    // size
    int textureH, textureW, texture_h, texture_w;
    int lossX{-1}, lossY{-1}, loss_x{-1}, loss_y{-1};

    int dpi_count{0};
    int res_h{-1}, res_w{-1};
    int *dpi{nullptr};
    int *dpi_n{nullptr};

    // const param
    const double min_overlap = MIN_DELTA;
    double half_min_delta = HALF_MIN_DELTA;
    double min_delta = SCALE_SIZE;
    double consume{0};

    double *pl_seam{nullptr};
    double *po_seam{nullptr};
    std::pair<int, int> *o_dpi{nullptr};

    // model choose
    void modelR();
    void modelGB();
    void modelLB();

    // generate
    void generateImg(int time_it);
    static int generateRand(int begin, int end);
    double generateRandD(double begin, double end);

    // check border situation
    static double len_v(const cv::Vec3i& v);
    static double len_V(const cv::Vec3i& vec);
    static bool edgeJudge(int he, int wi, int i, int j);

    // count num
    void countSeams();
    void countDpi() const;
    void countOldCut(int &OCx, int &OCy, int &OCX, int &OCY, int OCd);
    int countDpis(int r, int g, int b, int a);
    double countSeam(int i, int j, int edgeW, int edgeH);
    static double countLoss(const cv::Vec3i &c_ch1, const cv::Vec3i &c_ch2, const cv::Vec3i &C_ch1, const cv::Vec3i &C_ch2);

    // fft accelerate
    static cv::Mat fftAcc(Graph *ori, Graph* now);
    double dipRes();

    void reset(int h, int w);
    static int p2i(int t_w, int pi, int pj);

    MaxFlow *max_flow;
    Graph *source_img, *target_img, *out_img;

    // delta in x, y direction
    constexpr static const int delta_x[4] = {1, -1, 0, 0};
    constexpr static const int delta_y[4] = {0, 0, 1, -1};
};

GraphCutTexture::GraphCutTexture(std::string path) {
    max_flow = nullptr;

    target_img = nullptr;
    out_img =  nullptr;

    // read in
    source_img = new Graph(path);

    textureH = source_img->get_h();
    textureW = source_img->get_w();

    texture_h = textureH / 3;
    texture_w = textureW / 3;
}

GraphCutTexture::~GraphCutTexture() = default;

bool GraphCutTexture::graphCut(int h, int w, int i_time) {
    int iteration = i_time;
    int const_i = iteration;

    if (textureW > w || textureH > h || iteration < 0) {
        return false;
    }

    reset(h, w);

    int curr = 0;
    for (; dpi_count != res_w * res_h; ){
        countDpi();

        switch (model) {
            case RANDOM:
                modelR();
                break;
            case GLOBAL_BEST:
                modelGB();
                break;
            case LOCAL_BEST:
                modelLB();
                break;
        }

        curr++;
        generateImg(curr);

        std::cout << "\r";
        printf("Graph cut process: %.2lf%%", dpi_count * 100.0 / (res_h * res_w));
        std::cout << std::flush;
    }

    countDpi();
    printf("\n\n");

    while(iteration--){
        std::cout << "\r";
        printf("Refine process: %.2lf%%", (const_i - iteration) * 100.0 / const_i);
        std::cout << std::flush;

        switch (model) {
            case RANDOM:
                modelR();
                break;
            case GLOBAL_BEST:
                modelGB();
                break;
            case LOCAL_BEST:
                modelLB();
                break;
        }

        curr++;
        generateImg(curr);
    }
    printf("\n\n");

    return true;
}

cv::Mat GraphCutTexture::fftAcc(Graph *ori, Graph *now) {
    cv::Mat c_ch[3];
    cv::Mat C_ch[3];
    cv::Mat cha_c[3];
    cv::Mat cha_C[3];

    cv::split(ori->img, c_ch);
    cv::split(now->img, C_ch);

    for (int i = 0; i < 3; i++){
        cha_c[i] = cv::Mat_<float>(c_ch[i]);
        cha_C[i] = cv::Mat_<float>(C_ch[i]);
    }

    int d_w = cv::getOptimalDFTSize(cha_c[0].rows + cha_C[0].rows - 1);
    int d_h = cv::getOptimalDFTSize(cha_c[0].cols + cha_C[0].cols - 1);

    cv::Mat im_c = cv::Mat::zeros(d_w, d_h, CV_32F);

    for (int i = 0; i < 3; i++){
        cv::Mat d_af = cv::Mat::zeros(d_w, d_h, CV_32F);
        cv::Mat d_bf = cv::Mat::zeros(d_w, d_h, CV_32F);

        cv::Mat d_ap = d_af(cv::Rect(0, 0, cha_c[i].cols, cha_c[i].rows));
        cv::Mat d_bp = d_bf(cv::Rect(0, 0, cha_C[i].cols, cha_C[i].rows));

        cha_C[i].copyTo(d_bp);
        cha_c[i].copyTo(d_ap);

        cv::dft(d_af, d_af, 0, cha_c[0].rows);
        cv::dft(d_bf, d_bf, 0, cha_C[0].rows);

        cv::mulSpectrums(d_af, d_bf, d_af, 0, false);
        cv::idft(d_af, d_af, CV_HAL_DFT_SCALE, cha_c[i].rows + cha_C[i].rows - 1);

        im_c += d_af;
    }

    return im_c;
}

void GraphCutTexture::reset(int h, int w) {
    po_seam = new double[w * h]{0.0};
    pl_seam = new double[w * h]{0.0};

    dpi = new int[w * h];
    dpi_n = new int[w * h]{0};

    res_h = h;
    res_w = w;

    max_flow = new MaxFlow(3 * w * h, 8 * w * h);
    out_img = new Graph(h, w, source_img->img.type());

    o_dpi = new std::pair<int, int>[w * h];
    dpi_count = textureH * textureW;

    // init
    for (int i = 0; i < w * h; i++) {
        o_dpi[i] = std::make_pair(-1, -1);
        dpi[i] = -1;
    }

    for(int i = 0; i < textureH; i++){
        for(int j = 0; j < textureW; j++){
            int idx = p2i(w, i, j);

            o_dpi[idx] = std::make_pair(i, j);
            dpi[idx] = 0;

            out_img->img.at<cv::Vec3b>(i, j) = source_img->img.at<cv::Vec3b>(i, j);
        }
    }
}

// get img
void GraphCutTexture::generateImg(int time_it) {
    max_flow->reset();
    max_flow->set_begin(0);
    max_flow->set_end(1);

    int c_node = 1;
    int co_imh = target_img->get_h();
    int co_imw = target_img->get_w();

    std::map<int, int> i2p;
    std::vector<std::pair<int, int>> cover_p;

    for (int i = lossX; i < std::min(lossX + co_imh, res_h); i++){
        for (int j = lossY; j < std::min(lossY + co_imw, res_w); j++){
            int id = p2i(res_w, i, j);

            if (dpi[id] != -1){
                i2p[id] = ++c_node;
                cover_p.emplace_back(i, j);
            }
        }
    }

    double cons = consume;

    for (auto pos: cover_p){
        int index = p2i(res_w, pos.first, pos.second);
        int idx = i2p[index];

        int oub = 0;
        int teb = 0;

        for (int k = 0; k < 4; k++){
            int tw = pos.first + delta_x[k];
            int th = pos.second + delta_y[k];

            if (!edgeJudge(res_h, res_w, tw, th)){
                if (dpi_count != res_w * res_h && edgeJudge(co_imh, co_imw, tw - lossX, th - lossY)) {
                    teb++;
                }
                continue;
            }

            int n_idx = p2i(res_w, tw, th);
            if (i2p.count(n_idx)) {
                if (k & 1) {
                    continue;
                }
                cv::Vec3i patch_color_u = target_img->img
                        .at<cv::Vec3b>(pos.first - lossX, pos.second - lossY);
                cv::Vec3i patch_color_v = target_img->img
                        .at<cv::Vec3b>(tw - lossX, th - lossY);
                cv::Vec3i color_u_from_u = out_img->img.at<cv::Vec3b>(pos.first, pos.second);
                cv::Vec3i color_v_from_v = out_img->img.at<cv::Vec3b>(tw, th);

                if (dpi[index] == dpi[n_idx]) {
                    double w_cut = countLoss(color_u_from_u, color_v_from_v, patch_color_u, patch_color_v);

                    max_flow->insert(idx, i2p[n_idx], w_cut);
                }
                else {
                    c_node++;
                    cv::Vec3i col_r = source_img->img.at<cv::Vec3b>(o_dpi[index].first + delta_x[k], o_dpi[index].second + delta_y[k]);
                    cv::Vec3i col_g =source_img->img.at<cv::Vec3b>(o_dpi[n_idx].first - delta_x[k], o_dpi[n_idx].second - delta_y[k]);

                    double c_x = countLoss(color_u_from_u, col_r, patch_color_u, patch_color_v);
                    double c_X = countLoss(patch_color_u, patch_color_v, col_g, color_v_from_v);
                    double cXX = countLoss(color_u_from_u, col_r, col_g, color_v_from_v);

                    max_flow->insert(c_node, 1, cXX);
                    max_flow->insert(idx, c_node, c_x);
                    max_flow->insert(c_node, i2p[n_idx], c_X);

                    consume -= cXX;
                }
            }
            else {
                if (edgeJudge(co_imh, co_imw, tw - lossX, th - lossY)) {
                    teb++;
                }
                else if (dpi[n_idx] != -1) {
                    oub++;
                }
            }
        }
        if (teb + oub){
            if (teb >= oub) {
                max_flow->insert(idx, 1, INF);
            }
            else {
                max_flow->insert(0, idx, INF);
            }
        }
    }

    max_flow->set_n_num(c_node + 10);
    double max_flow_ans = max_flow->dinitz();

    if (max_flow_ans  >= INF) {
        consume = cons;
        return;
    }

    consume += max_flow_ans;

    bool *s_arrive = new bool[c_node + 10];
    max_flow->search(s_arrive);

    for (int i = lossX; i < std::min(lossX + co_imh, res_h); i++) {
        for (int j = lossY; j < std::min(lossY + co_imw, res_w); j++) {
            int idx = p2i(res_w, i, j);
            if ((dpi[idx] != -1 && !s_arrive[i2p[idx]]) || dpi[idx] == -1) {
                if (dpi[idx] == -1) {
                    dpi_count++;
                }

                out_img->img.at<cv::Vec3b>(i, j) = target_img->img.at<cv::Vec3b>(i - lossX, j - lossY);
                dpi[idx] = time_it;
                o_dpi[idx] = std::make_pair(i - lossX + loss_x, j - lossY + loss_y);
            }
        }
    }
}

// random int
int GraphCutTexture::generateRand(int begin, int end) {
    static std::default_random_engine e(time(nullptr));
    std::uniform_int_distribution<int> u(begin, end);

    return u(e);
}

// random double
double GraphCutTexture::generateRandD(double begin, double end) {
    static std::default_random_engine e(time(nullptr));
    std::uniform_real_distribution<double> u(begin, end);

    return u(e);
}

int GraphCutTexture::p2i(int t_w, int pi, int pj) {
    return pj + pi * t_w;
}

void GraphCutTexture::setModel(int model) {
    this->model = model;
}

void GraphCutTexture::countDpi() const {
    for (int i = 0; i < res_h; i++) {
        for (int j = 0; j < res_w; j++) {
            int idx = p2i(res_w, i, j);

            if (dpi[idx] == -1) {
                dpi_n[idx] = 0;
            }
            else {
                dpi_n[idx] = 1;
            }

            if (i != 0) {
                dpi_n[idx] += dpi_n[p2i(res_w, i - 1, j)];
            }

            if (j != 0) {
                dpi_n[idx] += dpi_n[p2i(res_w, i, j - 1)];
            }

            if (!(i == 0 || j == 0)) {
                dpi_n[idx] -= dpi_n[p2i(res_w, i - 1, j - 1)];
            }
        }
    }
}

int GraphCutTexture::countDpis(int r, int g, int b, int a) {
    int ret = dpi_n[p2i(res_w, b, a)] + ((r == 0 || g == 0) ? 0 : dpi_n[p2i(res_w, r - 1, g - 1)]);
    int tmp = (r == 0 ? 0 : dpi_n[p2i(res_w, r - 1, a)]) + (g == 0 ? 0 : dpi_n[p2i(res_w, b, g - 1)]);
    return ret - tmp;
}

double GraphCutTexture::dipRes() {
    int dpi_tc = 0;
    double arg = 0;

    for (int i = 0; i < res_h; i++){
        for (int j = 0; j < res_w; j++){
            int index = p2i(res_w, i, j);
            if (dpi[index] == -1) {
                continue;
            }

            dpi_tc++;
            arg += len_v(out_img->img.at<cv::Vec3b>(i, j));
        }
    }
    arg /= dpi_tc;

    double ret = 0;
    for (int i = 0; i < res_h; i++){
        for (int j = 0; j < res_w; j++){
            int index = p2i(res_w, i, j);
            if (dpi[index] == -1) {
                continue;
            }
            double ou_im = len_v(out_img->img.at<cv::Vec3b>(i, j)) - arg;

            ret += (ou_im * ou_im);
        }
    }

    ret /= dpi_tc;

    return ret;
}

void GraphCutTexture::modelR() {
    std::vector<std::pair<int, int>> img_random;

    for(int i = 0; i < res_h; i++){
        for(int j = 0; j < res_w; j++){
            int edge_w = std::min(i + textureH, res_h);
            int edge_h = std::min(j + textureW, res_w);

            int count_dpi_part = (edge_w - i) * (edge_h - j);
            int count_dpi_all = countDpis(i, j, edge_w - 1, edge_h - 1);

            if ((count_dpi_all < count_dpi_part * min_overlap) ||
                (count_dpi_part == count_dpi_all && dpi_count != res_w * res_h)) {
                continue;
            }

            img_random.emplace_back(i, j);
        }
    }

    loss_x = 0;
    loss_y = 0;

    int rand = generateRand(0, (int) (img_random.size()));
    lossX = img_random[rand].first;
    lossY = img_random[rand].second;

    target_img = new Graph(textureH, textureW, source_img->img.type());
    target_img->img = source_img->img.clone();
}

void GraphCutTexture::modelGB() {
    int check_loss_w, check_loss_h, check_loss_W, check_loss_H;
    countOldCut(check_loss_w, check_loss_h, check_loss_W, check_loss_H, 1);

    std::vector<std::pair<int, int>> vec_GB;
    for (int i = 0; i < check_loss_w + 1; i++){
        for(int j = 0; j < check_loss_h + 1; j++){
            int edge_w = std::min(i + textureH, res_h);
            int edge_h = std::min(j + textureW, res_w);

            if ((edge_w < check_loss_W) || (edge_h < check_loss_H)) {
                continue;
            }

            vec_GB.emplace_back(i, j);
        }
    }

    Graph *fft_acc_img = new Graph(res_h, res_w, 21);

    for (int i = 0; i < res_h; i++) {
        for(int j = 0; j < res_w; j++){
            if (i < textureH && j < textureW) {
                fft_acc_img->img.at<cv::Vec3f>(i, j) = source_img->img.at<cv::Vec3b>(i, j);
            }
            else {
                fft_acc_img->img.at<cv::Vec3f>(i, j) = cv::Vec3f(0, 0, 0);
            }
        }
    }

    cv::flip(fft_acc_img->img, fft_acc_img->img, -1);

    cv::Mat conv_ans = fftAcc(fft_acc_img, out_img);
    Graph *mask = new Graph(res_h, res_w, source_img->img.type());

    for (int i = 0; i < res_h; i++){
        for(int j = 0; j < res_w; j++){
            int id = p2i(res_w, i, j);
            if (dpi[id] == -1) {
                mask->img.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, 0);
            }
            else{
                mask->img.at<cv::Vec3b>(i, j) = cv::Vec3b(1, 1, 1);
            }

            cv::Vec3f acc_v = fft_acc_img->img.at<cv::Vec3f>(i, j);

            fft_acc_img->img.at<cv::Vec3f>(i, j) = cv::Vec3f(acc_v[0] * acc_v[0], acc_v[1] * acc_v[1], acc_v[2] * acc_v[2]);
        }
    }

    cv::Mat conv_ans2 = fftAcc(fft_acc_img, mask);

    double *c_count = new double[res_h * res_w];

    for(int i = 0; i < res_h; i++){
        for(int j = 0; j < res_w; j++){
            int idx = p2i(res_w, i, j);

            cv::Vec3i tout = out_img->img.at<cv::Vec3b>(i, j);

            c_count[idx] = (tout[0] * tout[0] + tout[1] * tout[1] + tout[2] * tout[2]);
            c_count[idx] += (i == 0 ? 0 : c_count[p2i(res_w, i - 1, j)]);
            c_count[idx] += (j == 0 ? 0 : c_count[p2i(res_w, i, j - 1)]);
            c_count[idx] -= ((i == 0 || j == 0) ? 0 : c_count[p2i(res_w, i - 1, j - 1)]);
        }
    }

    double var = dipRes();
    double al_num = 0;
    std::vector<double>count_node;

    for (auto pos: vec_GB){
        int pos_w = pos.first;
        int pos_h = pos.second;

        int edge_w = std::min(res_h, pos_w + textureH);
        int edge_h = std::min(res_w, pos_h + textureW);

        double tmp = (pos_w == 0 ? 0 : c_count[p2i(res_w, pos_w - 1, edge_h - 1)]) + (pos_h == 0 ? 0 : c_count[p2i(res_w, edge_w - 1, pos_h - 1)]);
        double c_count_O = c_count[p2i(res_w, edge_w - 1, edge_h - 1)] + ((pos_w == 0 || pos_h == 0) ? 0 : c_count[p2i(res_w, pos_w - 1, pos_h - 1)]);
        c_count_O -= tmp;
        double c_count_I = conv_ans2.at<float>(res_h - 1 + pos_w, res_w - 1 + pos_h);
        double c_count_T = 2 * conv_ans.at<float>(res_h - 1 + pos_w, res_w - 1 + pos_h);

        double offset = c_count_O - c_count_T + c_count_I;
        int cover = countDpis(pos_w, pos_h, edge_w - 1, edge_h - 1);
        offset /= cover;

        double al = min_delta * exp(-offset / (half_min_delta * var));
        al_num += al;

        count_node.push_back(al);
    }
    double al_cho = generateRandD(0, al_num);

    int count_p_n = 0;

    lossX = -1;
    lossY = -1;

    for (auto pos: vec_GB){
        al_cho -= count_node[count_p_n];
        if (al_cho <= 0) {
            lossX = pos.first;
            lossY = pos.second;

            break;
        }
        count_p_n++;
    }

    loss_x = 0;
    loss_y = 0;

    target_img = new Graph(textureH, textureW, source_img->img.type());
    target_img->img = source_img->img.clone();
}

void GraphCutTexture::modelLB() {
    int check_loss_w, check_loss_h, check_loss_W, check_loss_H;
    countOldCut(check_loss_w, check_loss_h,check_loss_W, check_loss_H, 3);

    int t_H = texture_h;
    int t_W = texture_w;

    if (check_loss_w > 0){
        t_H++;
        check_loss_w--;
    }

    if(check_loss_h > 0){
        t_W++;
        check_loss_h--;
    }

    if (check_loss_w + t_H < res_h) {
        t_H++;
    }

    if (check_loss_h + t_W < res_w) {
        t_W++;
    }

    double c_i = 0;
    int cover_c = 0;
    Graph *fft_acc_g = new Graph(textureH, textureW, 21);
    Graph *cove = new Graph(textureH, textureW, 21);

    for (int i = 0; i < textureH; i++){
        for(int j = 0; j < textureW; j++){
            if(i < t_H && j < t_W){
                int nw = i + check_loss_w;
                int nh = j + check_loss_h;
                int nxt_id = p2i(res_w, nw, nh);

                if (dpi[nxt_id] == -1) {
                    continue;
                }

                cover_c++;

                cv::Vec3f o_img = out_img->img.at<cv::Vec3b>(nw, nh);
                c_i += (o_img[0] * o_img[0]) + (o_img[1] * o_img[1]) + (o_img[2] * o_img[2]);

                fft_acc_g->img.at<cv::Vec3f>(i, j) = o_img;
                cove->img.at<cv::Vec3f>(i, j) = cv::Vec3f(1, 1, 1);
            }
        }
    }
    cv::flip(fft_acc_g->img, fft_acc_g->img, -1);
    cv::flip(cove->img, cove->img, -1);

    cv::Mat conv_ans = fftAcc(fft_acc_g, source_img);

    for (int i = 0; i < textureH; i++) {
        for (int j = 0; j < textureW; j++) {
            cv::Vec3f src_img = source_img->img.at<cv::Vec3b>(i, j);
            fft_acc_g->img.at<cv::Vec3f>(i, j) = cv::Vec3f(src_img[0] * src_img[0], src_img[1] * src_img[1], src_img[2] * src_img[2]);
        }
    }

    cv::Mat ftt_acc_res = fftAcc(fft_acc_g, cove);

    int w_h = textureH - t_H;
    int h_h = textureW - t_W;

    loss_x = -1;
    loss_y = -1;

    double res_dpi = dipRes();
    double delta_o_sum = 0;
    std::vector<double>delta_o_record;

    for (int i = 0; i < w_h; i++){
        for (int j = 0; j < h_h; j++) {
            double c_ftt_c = ftt_acc_res.at<float>(textureH - 1 + i, textureW - 1 + j);
            double ftt_o_c = 2.0 * conv_ans.at<float>(textureH - 1 + i, textureW - 1 + j);

            double offset = c_i - ftt_o_c + c_ftt_c;
            offset /= cover_c;

            double delta_o = min_delta * exp(-offset / (half_min_delta * res_dpi));

            delta_o_sum += delta_o;

            delta_o_record.push_back(delta_o);
        }
    }
    double delta_cho = generateRandD(0, delta_o_sum);

    bool che_cho = false;
    int delta_o_c = 0;

    for (int i = 0; i < w_h; i++) {
        for (int j = 0; j < h_h; j++){
            delta_cho -= delta_o_record[delta_o_c];

            if (delta_cho <= 0) {
                che_cho = true;

                loss_x = i;
                loss_y = j;

                break;
            }
            delta_o_c++;
        }
        if (che_cho) {
            break;
        }
    }

    lossX = check_loss_w;
    lossY = check_loss_h;

    target_img = new Graph(t_H, t_W, source_img->img.type());
    target_img->img = source_img->img(cv::Rect(loss_y, loss_x, t_W, t_H)).clone();
}

double GraphCutTexture::countLoss(const cv::Vec3i &c_ch1, const cv::Vec3i &c_ch2, const cv::Vec3i &C_ch1, const cv::Vec3i &C_ch2) {
    return len_v(c_ch1 - C_ch1) + len_v(c_ch2 - C_ch2);
}

double GraphCutTexture::len_v(const cv::Vec3i& vec) {
    return sqrt(len_V(vec));
}

double GraphCutTexture::len_V(const cv::Vec3i& vec) {
    return (double)(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

bool GraphCutTexture::edgeJudge(int he, int wi, int i, int j) {
    return i >= 0 && i < he && j >= 0 && j < wi;
}

void GraphCutTexture::saveTexture(const std::string path) {
    cv::imwrite(path, out_img->img);
}

void GraphCutTexture::saveSeams(const std::string path_s) {
    Graph *seam = new Graph(res_h, res_w, out_img->img.type());

    seam->img = out_img->img.clone();

    for (int i = 0; i < res_h; i++){
        for (int j = 0; j < res_w; j++){
            bool che = false;
            int idx = p2i(res_w, i, j);

            for (int k = 0; k < 4; k++){
                int tw = delta_x[k] + i;
                int th = delta_y[k] + j;

                if (!edgeJudge(res_h, res_w, tw, th)) {
                    continue;
                }

                int t_idx = p2i(res_w, tw, th);
                if (dpi[t_idx] != dpi[idx]) {
                    che = true;
                    break;
                }
            }

            if (che) {
                seam->img.at<cv::Vec3b>(i, j) = cv::Vec3b(205, 133, 63);
            }
        }
    }

    cv::imwrite(path_s, seam->img);
}

void GraphCutTexture::countSeams() {
    // count
    for (int i = 0; i < res_h; i++){
        for(int j = 0; j < res_w; j++){
            int idx = p2i(res_w, i, j);
            pl_seam[idx] = po_seam[idx] = 0;

            if (dpi[idx] == -1) {
                continue;
            }

            for (int k = 1; k < 4; k += 2){
                int tw = i + delta_x[k];
                int th = j + delta_y[k];

                if (!edgeJudge(res_h, res_w, tw, th)) {
                    continue;
                }

                int t_idx = p2i(res_w, tw, th);
                if (dpi[t_idx] == -1) {
                    continue;
                }

                if (dpi[t_idx] == dpi[idx]) {
                    continue;
                }

                cv::Vec3i c_r = out_img->img.at<cv::Vec3b>(i, j);
                cv::Vec3i c_g = out_img->img.at<cv::Vec3b>(tw, th);
                cv::Vec3i c_b = source_img->img.at<cv::Vec3b>(o_dpi[idx].first + delta_x[k], o_dpi[idx].second + delta_y[k]);
                cv::Vec3i c_a = source_img->img.at<cv::Vec3b>(o_dpi[t_idx].first - delta_x[k], o_dpi[t_idx].second - delta_y[k]);

                double w_cut = countLoss(c_r, c_b, c_a, c_g);

                if (k == 1) {
                    pl_seam[idx] = w_cut;
                }
                else {
                    po_seam[idx] = w_cut;
                }
            }
        }
    }

    for (int i = 0; i < res_h; i++){
        for(int j = 0; j < res_w; j++){
            int index = p2i(res_w, i, j);

            pl_seam[index] += (i == 0 ? 0 : pl_seam[p2i(res_w, i - 1, j)]);
            pl_seam[index] += (j == 0 ? 0 : pl_seam[p2i(res_w, i, j - 1)]);
            pl_seam[index] -= ((i == 0 || j == 0) ? 0 : pl_seam[p2i(res_w, i - 1, j - 1)]);
            po_seam[index] += (i == 0 ? 0 : po_seam[p2i(res_w, i - 1, j)]);
            po_seam[index] += (j == 0 ? 0 : po_seam[p2i(res_w, i, j - 1)]);
            po_seam[index] -= ((i == 0 || j == 0) ? 0 : po_seam[p2i(res_w, i - 1, j - 1)]);
        }
    }
}

double GraphCutTexture::countSeam(int i, int j, int edgeW, int edgeH) {
    double tmp_l = pl_seam[p2i(res_w, i, edgeH)] + (j == 0? 0: pl_seam[p2i(res_w, edgeW, j - 1)]);
    double tmp_L = pl_seam[p2i(res_w, edgeW, edgeH)] + (j == 0 ? 0: pl_seam[p2i(res_w, i, j - 1)]);

    double tmp_o = po_seam[p2i(res_w, edgeW, j)] + (i == 0? 0: po_seam[p2i(res_w, i - 1, edgeH)]);
    double tmp_O = po_seam[p2i(res_w, edgeW, edgeH)] + (i == 0? 0: po_seam[p2i(res_w, i - 1, j)]);

    return tmp_L + tmp_O - tmp_l - tmp_o;
}

void GraphCutTexture::countOldCut(int &OCx, int &OCy, int &OCX, int &OCY, int OCd) {
    countSeams();

    double check_loss = -1;
    int add_dpi = 0;

    int uw = res_h - texture_h / 2;
    int uh = res_w - texture_w / 2;

    for (int i = 0; i < uw; i++) {
        for (int j = 0; j < uh; j++) {
            int edgeW = std::min(i + texture_h, res_h);
            int edgeH = std::min(j + texture_w, res_w);

            int dpi_s = countDpis(i, j, edgeW - 1, edgeH - 1);
            int dpi_c = (edgeW - i) * (edgeH - j);

            if (dpi_s < OCd * dpi_c * min_overlap) {
                continue;
            }

            if (dpi_count != res_w * res_h && dpi_c == dpi_s) {
                continue;
            }

            double del_s = countSeam(i, j, edgeW - 1, edgeH - 1);
            int a_dpi = dpi_c - dpi_s;

            if (del_s > check_loss || (del_s == check_loss && a_dpi > add_dpi)) {
                OCX = edgeW;
                OCY = edgeH;

                check_loss = del_s;
                add_dpi = a_dpi;

                OCx = i;
                OCy = j;
            }
        }
    }
}

int GraphCutTexture::getModel() {
    return this->model;
}

int GraphCutTexture::getWidth() {
    return this->textureW;
}

int GraphCutTexture::getHeight() {
    return this->textureH;
}

#endif //GRAPHCUT_GRAPHCUTTEXTURE_HPP
