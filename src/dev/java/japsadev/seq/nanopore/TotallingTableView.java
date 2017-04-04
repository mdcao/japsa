package japsadev.seq.nanopore;

import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javafx.application.Application;
import javafx.beans.Observable;
import javafx.beans.binding.Bindings;
import javafx.beans.binding.IntegerBinding;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;
import javafx.util.StringConverter;
import javafx.util.converter.DefaultStringConverter;
import javafx.util.converter.NumberStringConverter;

public class TotallingTableView extends Application {

    @Override
    public void start(Stage primaryStage) {

        TableView<Item> table = new TableView<>();
        table.setEditable(true);

        table.getColumns().add(
                column("Item", Item::nameProperty, new DefaultStringConverter()));

        table.getColumns().add(
                column("Value", Item::valueProperty, new NumberStringConverter()));

        table.setItems(FXCollections.observableArrayList(
                item -> new Observable[] {item.valueProperty() }));

        IntStream.rangeClosed(1, 20)
            .mapToObj(i -> new Item("Item "+i, i))
            .forEach(table.getItems()::add);

        IntegerBinding total = Bindings.createIntegerBinding(() -> 
            table.getItems().stream().collect(Collectors.summingInt(Item::getValue)),
            table.getItems());

        Label totalLabel = new Label();
        totalLabel.textProperty().bind(Bindings.format("Total: %d", total));

        Button add = new Button("Add item");
        add.setOnAction(e -> 
            table.getItems().add(new Item("New Item", table.getItems().size() + 1)));

        Button remove = new Button("Remove");
        remove.disableProperty().bind(
                Bindings.isEmpty(table.getSelectionModel().getSelectedItems()));

        remove.setOnAction(e -> 
            table.getItems().remove(table.getSelectionModel().getSelectedItem()));

        HBox buttons = new HBox(5, add, remove);
        buttons.setAlignment(Pos.CENTER);
        VBox controls = new VBox(5, totalLabel, buttons);
        VBox.setVgrow(totalLabel, Priority.ALWAYS);
        totalLabel.setMaxWidth(Double.MAX_VALUE);
        totalLabel.setAlignment(Pos.CENTER_RIGHT);

        BorderPane root = new BorderPane(table, null, null, controls, null);
        Scene scene = new Scene(root, 800, 600);
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    private <S,T> TableColumn<S,T> column(String title, 
            Function<S, ObservableValue<T>> property, StringConverter<T> converter) {
        TableColumn<S,T> col = new TableColumn<>(title);
        col.setCellValueFactory(cellData -> property.apply(cellData.getValue()));

        col.setCellFactory(TextFieldTableCell.forTableColumn(converter));

        return col ;
    }

    public static class Item {
        private final StringProperty name = new SimpleStringProperty();
        private final IntegerProperty value = new SimpleIntegerProperty();

        public Item(String name, int value) {
            setName(name);
            setValue(value);
        }

        public final StringProperty nameProperty() {
            return this.name;
        }

        public final java.lang.String getName() {
            return this.nameProperty().get();
        }

        public final void setName(final java.lang.String name) {
            this.nameProperty().set(name);
        }

        public final IntegerProperty valueProperty() {
            return this.value;
        }

        public final int getValue() {
            return this.valueProperty().get();
        }

        public final void setValue(final int value) {
            this.valueProperty().set(value);
        }


    }

    public static void main(String[] args) {
        launch(args);
    }
}
