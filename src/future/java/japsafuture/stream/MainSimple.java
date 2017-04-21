package japsafuture.stream;

import java.util.Arrays;
import java.util.concurrent.SubmissionPublisher;

/**
 * Created by minhduc on 21/04/17.
 */
public class MainSimple {

    public static void main(String[] args) throws InterruptedException {


        //Create Publisher
        SubmissionPublisher<String> publisher = new SubmissionPublisher<>();

        //Register Subscriber
        MySubscriberSimple<String> subscriber = new MySubscriberSimple<>();
        publisher.subscribe(subscriber);

        //Publish items
        System.out.println("Publishing Items...");
        String[] items = {"1", "x", "2", "x", "3", "x"};
        Arrays.asList(items).stream().forEach(i -> publisher.submit(i));
        publisher.close();

    }
}