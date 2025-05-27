package knu.lsy.shapeoverlap.exception;

import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.ControllerAdvice;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.validation.BindException;
import jakarta.validation.ConstraintViolationException;

import java.util.HashMap;
import java.util.Map;

@ControllerAdvice
public class GlobalExceptionHandler {

    @ExceptionHandler(ConstraintViolationException.class)
    public ResponseEntity<Map<String, Object>> handleValidationExceptions(ConstraintViolationException ex) {
        Map<String, Object> response = new HashMap<>();
        response.put("success", false);
        response.put("error", "Validation error");
        response.put("details", ex.getMessage());
        response.put("timestamp", System.currentTimeMillis());
        
        return ResponseEntity.badRequest().body(response);
    }

    @ExceptionHandler(BindException.class)
    public ResponseEntity<Map<String, Object>> handleBindExceptions(BindException ex) {
        Map<String, Object> response = new HashMap<>();
        response.put("success", false);
        response.put("error", "Invalid request parameters");
        response.put("details", ex.getBindingResult().getAllErrors());
        response.put("timestamp", System.currentTimeMillis());
        
        return ResponseEntity.badRequest().body(response);
    }

    @ExceptionHandler(Exception.class)
    public ResponseEntity<Map<String, Object>> handleAllExceptions(Exception ex) {
        Map<String, Object> response = new HashMap<>();
        response.put("success", false);
        response.put("error", "Internal server error");
        response.put("message", ex.getMessage());
        response.put("timestamp", System.currentTimeMillis());
        
        return ResponseEntity
            .status(HttpStatus.INTERNAL_SERVER_ERROR)
            .body(response);
    }
} 