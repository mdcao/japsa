package japsafuture.stream;

import static java.lang.Thread.currentThread;

import java.util.Random;
import java.util.concurrent.Flow.Subscriber;
import java.util.concurrent.Flow.Subscription;

public class MySubscriber implements Subscriber<Integer> {

    private static final String LOG_MESSAGE_FORMAT = "Subscriber %s >> [%s] %s%n";

    private static final int DEMAND = 3;
    private static final Random RANDOM = new Random();

    private String name;
    private Subscription subscription;

    private int count;

    public MySubscriber(String name) {
        this.name = name;
    }

    @Override
    public void onSubscribe(Subscription subscription) {
        log("Subscribed");
        this.subscription = subscription;

        count = DEMAND;
        requestItems(DEMAND);
    }

    private void requestItems(int n) {
        log("Requesting %d new items...", n);
        subscription.request(n);
    }

    @Override
    public void onNext(Integer item) {
        if (item != null) {
            log(item.toString());

            synchronized (this) {
                count--;

                if (count == 0) {
                    if (RANDOM.nextBoolean()) {
                        count = DEMAND;
                        requestItems(count);
                    } else {
                        count = 0;
                        log("Cancelling subscription...");
                        subscription.cancel();
                    }
                }
            }
        } else {
            log("Null Item!");
        }
    }

    @Override
    public void onComplete() {
        log("Complete!");
    }

    @Override
    public void onError(Throwable t) {
        log("Subscriber Error >> %s", t);
    }

    private void log(String message, Object... args) {
        String fullMessage = String.format(LOG_MESSAGE_FORMAT, this.name, currentThread().getName(), message);

        System.out.printf(fullMessage, args);
    }
}
