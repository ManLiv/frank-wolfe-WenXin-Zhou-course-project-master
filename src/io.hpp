#ifndef IO_HPP_
#define IO_HPP_

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

#include <fstream>

typedef enum {
    UNKNOWN_METADATA, NUMBER_OF_ZONES, NUMBER_OF_NODES, FIRST_THRU_NODE, NUMBER_OF_LINKS, TOTAL_OD_FLOW, LOCATION, END_OF_METADATA, NUMBER_OF_TOLLS
} meta_data_label;

std::pair<meta_data_label, std::string> read_meta_data(const std::string& str) {

    std::vector<std::string> result; // #2: Search for tokens

    boost::split(result, str, boost::is_any_of("<>"), boost::token_compress_on);
    boost::trim(result[2]);

    meta_data_label label(UNKNOWN_METADATA);

    if (result[1] == "NUMBER OF ZONES")
        label = NUMBER_OF_ZONES;
    if (result[1] == "NUMBER OF NODES")
        label = NUMBER_OF_NODES;
    if (result[1] == "NUMBER OF LINKS")
        label = NUMBER_OF_LINKS;
    if (result[1] == "FIRST THRU NODE")
        label = FIRST_THRU_NODE;
    if (result[1] == "TOTAL OD FLOW")
        label = TOTAL_OD_FLOW;
    if (result[1] == "END OF METADATA")
        label = END_OF_METADATA;
    return std::make_pair(label, result[2]);
}

template<typename graph_type, typename zone_list_type>
void load_network(const std::string& network_filename, graph_type &g, zone_list_type& centroids, int &num_centroids, bool & all_centroids) {

    typedef typename graph_type::vertex_descriptor v_type;

    int first_thru_node, num_arcs = 0, num_nodes = 0;
    std::ifstream network_file(network_filename.c_str());
    if (!network_file) {
        std::cout << "Network file does not exist!" << std::endl;
        exit(-1);
    }

    std::string line;
    bool loop = true;

    while (std::getline(network_file, line) && loop) {
        boost::trim(line);
        meta_data_label meta_data;
        std::string value;

        if (line[0] == '<') {
            boost::tie(meta_data, value) = read_meta_data(line);

            switch (meta_data) {
            case END_OF_METADATA:
                loop = false;
                break;
            case NUMBER_OF_ZONES:
                num_centroids = boost::lexical_cast<int>(value);
                break;
            case NUMBER_OF_NODES:
                num_nodes = boost::lexical_cast<int>(value);
                break;
            case NUMBER_OF_LINKS:
                num_arcs = boost::lexical_cast<int>(value);
                break;
            case FIRST_THRU_NODE:
                first_thru_node = boost::lexical_cast<int>(value);
                all_centroids = (first_thru_node == 1) ? true : false;
                break;
            default:
                break;
            }
        }
    }

    if (num_nodes > num_arcs) {
        std::cerr << "Fatal error! The graph is not connected!" << std::endl;
        exit(-1);
    }

    centroids.clear();
    for (int ii = 0; ii < num_nodes; ii++) {
        v_type u = boost::add_vertex(g);

        if (ii < num_centroids) {
            g[u].centroid = true;
            centroids.push_back(u);
        }
        else {
            g[u].centroid = false;
        }
    }

    
    while (std::getline(network_file, line)) {
        boost::trim(line);
        if (line[0] == '~' || line[0] == '<' || line.empty())
            continue;

        int source, destination, type, speed_limit, toll;
        double fft, B, length, capacity, power;

        std::stringstream ss(line);
        ss >> source >> destination >> capacity >> length >> fft >> B >> power >> speed_limit >> toll >> type;
        v_type source0 = source - 1;
        v_type destination0 = destination - 1;

        typename graph_type::edge_descriptor e = boost::add_edge(boost::vertex(source0, g), boost::vertex(destination0, g), g).first;

        g[e].cost_fun.initialize(capacity, fft, B, power, length, toll);
    }

    network_file.close();
}

template<typename matrix_type>
int load_trips(const std::string& trips_filename, matrix_type &D, std::vector<uint>& destination_count, double& totalDemand) {
    int num_centroids = 0, num_pairs = 0;
    double total_flow = 0.;

    std::ifstream trips_file(trips_filename.c_str());
    if (!trips_file) {
        std::cout << "Trips file does not exist!" << std::endl;
        exit(-1);
    }

    bool loop = true;
    std::string line;

    while (std::getline(trips_file, line) && loop) {
        boost::trim(line);
        if (line[0] == '~' || line.empty())
            continue;

        if (line[0] == '<') {
            meta_data_label meta_data;
            std::string value;
            boost::tie(meta_data, value) = read_meta_data(line);
            switch (meta_data) {
            case END_OF_METADATA:
                loop = false;
                break;
            case NUMBER_OF_ZONES:
                num_centroids = boost::lexical_cast<int>(value);
                break;
            case TOTAL_OD_FLOW:
                total_flow = boost::lexical_cast<double>(value);
                totalDemand = total_flow;
                break;
            default:
                break;
            }

        }
    }

    D = matrix_type(num_centroids, num_centroids, 0.);
    destination_count.resize(num_centroids, 0);

    int origin = -1;
    while (std::getline(trips_file, line)) {
        boost::trim(line);
        if (line.empty() || line[0] == '~' || line[0] == '<')
            continue;

        if (line[0] == 'O') {
            std::vector<std::string> result;
            boost::split(result, line, boost::is_any_of(" "), boost::token_compress_on);
            boost::trim(result[1]);

            // Origine
            origin = boost::lexical_cast<int>(result[1]) - 1;
            continue;
        }

        std::istringstream iss(line);
        std::string app;
        while (getline(iss, app, ';')) {
            boost::trim(app);
            if (app.empty())
                continue;

            std::vector<std::string> result;
            boost::split(result, app, boost::is_any_of(":"), boost::token_compress_on);

            boost::trim(result[0]);
            boost::trim(result[1]);
            int destination = boost::lexical_cast<int>(result[0]) - 1;
            double flow = boost::lexical_cast<double>(result[1]);

            // Non si considerano le domande di trasporto intra-zonali
            if (origin != destination && flow > 0.0) {
                D(origin, destination) = flow;
                destination_count[origin]++;
                num_pairs++;
            }
        }
    }

    trips_file.close();

    return num_pairs;
}


#endif /*IO_HPP_*/