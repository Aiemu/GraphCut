#include <iostream>
#include <GraphCutTexture.hpp>

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cout << "Usage: ./GraphCut <model> <iteration_time> <patch_path> <output_path> <output_width> <output_height>\n\n";
        std::cout << "e.g. ./GraphCut 0 100 green.gif green-out.png 128 128\n";
        return 0;
    }

    GraphCutTexture *GCT = new GraphCutTexture(argv[3]);
    int mo;
    sscanf(argv[1], "%d", &mo);
    GCT->setModel(mo);

    int w, h;
    sscanf(argv[5], "%d", &w);
    sscanf(argv[6], "%d", &h);

    int i_time;
    sscanf(argv[2], "%d", &i_time);
    if (GCT->graphCut(w, h, i_time)) {
        printf("Please check result at output path.\n");
    }

    std::string base_out_path("./output/");
    std::string output_path_ori("");
    std::string output_path_old_cut("seams_");
    std::string output_path(argv[4]);

    output_path_ori = base_out_path + output_path;
    GCT->saveTexture(output_path_ori);

    output_path_old_cut = base_out_path + output_path_old_cut + output_path;
    GCT->saveSeams(output_path_old_cut);

    return 0;
}
