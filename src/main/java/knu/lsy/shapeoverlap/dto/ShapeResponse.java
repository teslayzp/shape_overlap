package knu.lsy.shapeoverlap.dto;

public class ShapeResponse<T> {
    private String message;
    private T data;

    public ShapeResponse(String message, T data) {
        this.message = message;
        this.data = data;
    }

    // Getters and setters
    public String getMessage() { return message; }
    public void setMessage(String message) { this.message = message; }
    public T getData() { return data; }
    public void setData(T data) { this.data = data; }
}