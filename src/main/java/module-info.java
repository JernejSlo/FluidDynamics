module com.example.fluiddynamicsjavafx {
    requires javafx.controls;
    requires javafx.fxml;


    opens com.example.fluiddynamicsjavafx to javafx.fxml;
    exports com.example.fluiddynamicsjavafx;
}
