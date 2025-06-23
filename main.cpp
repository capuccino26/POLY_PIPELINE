#include <gtkmm.h>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <array>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <ctime>

class AppWindow : public Gtk::Window {
public:
    AppWindow() {
        set_title("STOMics Standalone Client");
        set_default_size(750, 550);

        button1.set_label("Test 1 (Python)");
        button2.set_label("Test 2 (Bash)");
        button3.set_label("Test 3 (R)");

        button1.signal_clicked().connect(sigc::mem_fun(*this, &AppWindow::on_analysis1));
        button2.signal_clicked().connect(sigc::mem_fun(*this, &AppWindow::on_analysis2));
        button3.signal_clicked().connect(sigc::mem_fun(*this, &AppWindow::on_analysis3));

        scroll.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
        scroll.add(output_view);
        output_view.set_editable(false);
        output_view.set_wrap_mode(Gtk::WRAP_WORD_CHAR);

        vbox.set_spacing(10);
        vbox.set_margin_top(10);
        vbox.set_margin_bottom(10);
        vbox.set_margin_left(10);
        vbox.set_margin_right(10);
        vbox.pack_start(button1, Gtk::PACK_SHRINK);
        vbox.pack_start(button2, Gtk::PACK_SHRINK);
        vbox.pack_start(button3, Gtk::PACK_SHRINK);
        vbox.pack_start(scroll, Gtk::PACK_EXPAND_WIDGET);

        add(vbox);
        show_all_children();
    }

protected:
    void on_analysis1() {
        std::string output = run_command("bash -c '~/anaconda3/bin/conda run -n poly_pipeline python scripts/script1.py'", "Python Analysis 1");
        output_view.get_buffer()->set_text(output);
    }

    void on_analysis2() {
        std::string output = run_command("bash scripts/script2.sh", "Bash Analysis 2");
        output_view.get_buffer()->set_text(output);
    }

    void on_analysis3() {
        std::string output = run_command("bash -c '~/anaconda3/bin/conda run -n poly_pipeline Rscript scripts/script3.R'", "R Analysis 3");
        output_view.get_buffer()->set_text(output);
    }

    std::string run_command(const std::string& cmd, const std::string& analysis_name) {
        std::array<char, 256> buffer;
        std::stringstream result;
        
        // Registro de tempo de execução
        auto start_time = std::chrono::system_clock::now();
        std::time_t start_c = std::chrono::system_clock::to_time_t(start_time);

        // Executa o comando
        FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
        if (!pipe) {
            result << "❌ Error executing: " << cmd << "\n";
            return result.str();
        }

        while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
            result << buffer.data();
        }
        pclose(pipe);

        // Calcula tempo de execução
        auto end_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::time_t end_c = std::chrono::system_clock::to_time_t(end_time);

        // Criação da entrada de log
        std::stringstream log_entry;
        log_entry << "===============================================\n";
        log_entry << "📅 Start: " << std::put_time(std::localtime(&start_c), "%Y-%m-%d %H:%M:%S") << "\n";
        log_entry << "⚙️  Command: " << cmd << "\n";
        log_entry << "📝 Analysis: " << analysis_name << "\n\n";
        log_entry << result.str() << "\n";
        log_entry << "⏱️  Duration: " << duration << "s\n";
        log_entry << "===============================================\n\n";

        // Salva o log no arquivo
        std::ofstream log_file("log.txt", std::ios::app);
        if (log_file.is_open()) {
            log_file << log_entry.str();
            log_file.close();
        }

        return result.str();
    }

private:
    Gtk::Box vbox{Gtk::ORIENTATION_VERTICAL};
    Gtk::Button button1, button2, button3;
    Gtk::ScrolledWindow scroll;
    Gtk::TextView output_view;
};

int main(int argc, char *argv[]) {
    auto app = Gtk::Application::create(argc, argv, "com.exemplo.analises");
    AppWindow window;
    return app->run(window);
}
