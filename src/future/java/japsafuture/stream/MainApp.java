package japsafuture.stream;

import java.util.Arrays;
import java.util.concurrent.SubmissionPublisher;

/**
 * Created by minhduc on 21/04/17.
 */
public class MainApp {

    public static void main(String[] args) throws InterruptedException {

        SubmissionPublisher<String> publisher1 = new SubmissionPublisher<>();

        /****************************************
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
         /****************************************/

        MyPublisher publisher = new MyPublisher();
        MySubscriber subscriberA = new MySubscriber("A");
        MySubscriber subscriberB = new MySubscriber("B");

        publisher.subscribe(subscriberA);
        publisher.subscribe(subscriberB);

        //publisher.waitUntilTerminated();

    }
}